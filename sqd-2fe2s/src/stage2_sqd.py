"""Stage 2 - sample-based quantum diagonalization.

Direct port of the IBM SQD chemistry tutorial
(https://quantum.cloud.ibm.com/docs/en/tutorials/sample-based-quantum-diagonalization):

    |HF>  ->  LUCJ (ffsim, CCSD t1/t2)  ->  sample  ->  configuration recovery
          ->  self-consistent subspace diagonalization

One deviation: the tutorial samples an LUCJ circuit on a QPU.  We sample it with
ffsim's exact sampler and use its ``global_depolarizing`` knob to stand in for
hardware noise, which is what makes configuration recovery do any work.  The
uniform-random-bitstring control is taken from the package quickstart.
"""

from __future__ import annotations

import argparse
import pathlib
import sys
import time
from functools import partial

import ffsim
import numpy as np
from pyscf import fci
from qiskit import QuantumCircuit, QuantumRegister
from qiskit.primitives.containers import BitArray
from qiskit_addon_sqd.counts import generate_bit_array_uniform
from qiskit_addon_sqd.fermion import (
    diagonalize_fermionic_hamiltonian,
    solve_sci,
    solve_sci_batch,
)

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import ffsim_patch  # noqa: E402,F401  (repairs to_parameters / optimize=True)
from common import RESULTS, load_json, save_json  # noqa: E402

# The tutorial uses max_cycle=200.  On this system that leaves the selected-CI Davidson
# short of convergence by 0.112 mHa even when diagonalizing the *entire* CAS space,
# which would masquerade as an SQD sampling error.  500 cycles reproduces CASCI to eight
# decimals; the "solver floor" row in the results reports whatever this setting costs.
SCI_MAX_CYCLE = 500


def build_circuit(t1, t2, norb, nelec, n_reps, optimize):
    """LUCJ ansatz, exactly as in the tutorial: nearest-neighbour same-spin
    interactions and all-to-all opposite-spin interactions."""
    pairs_aa = [(p, p + 1) for p in range(norb - 1)]
    pairs_ab = None
    ucj_op = ffsim.UCJOpSpinBalanced.from_t_amplitudes(
        t2=t2, t1=t1, n_reps=n_reps, interaction_pairs=(pairs_aa, pairs_ab),
        optimize=optimize,
    )
    qubits = QuantumRegister(2 * norb, name="q")
    circuit = QuantumCircuit(qubits)
    circuit.append(ffsim.qiskit.PrepareHartreeFockJW(norb, nelec), qubits)
    circuit.append(ffsim.qiskit.UCJOpSpinBalancedJW(ucj_op), qubits)
    circuit.measure_all()
    return circuit, ucj_op


def ansatz_energy(ucj_op, hcore, eri, e_nuc, norb, nelec):
    vec = ffsim.hartree_fock_state(norb, nelec)
    vec = ffsim.apply_unitary(vec, ucj_op, norb=norb, nelec=nelec)
    ham = ffsim.linear_operator(
        ffsim.MolecularHamiltonian(hcore, eri, constant=e_nuc), norb=norb, nelec=nelec
    )
    prob = np.abs(vec) ** 2
    prob /= prob.sum()
    return {
        "e_ansatz": float(np.real(np.vdot(vec, ham @ vec))),
        "participation_ratio": float(1.0 / np.sum(prob**2)),
        "max_probability": float(prob.max()),
    }


def sample_bitstrings(circuit, shots, norb, nelec, backend, depol, seed):
    """Sample the ansatz circuit on one of several simulator backends.

    ffsim      exact statevector sampling in the fixed-particle-number sector, with an
               optional global depolarizing channel.  Cheapest, and the only one that
               scales past ~30 qubits here, so it is the default.
    aer        qiskit-aer statevector with a gate-level depolarizing noise model
               (1-qubit `depol`, 2-qubit 10x `depol`) - closest stand-in for a device.
    aer-mps    same, matrix-product-state method; the option to reach for when the
               qubit count outgrows a dense statevector.

    Returns (BitArray, info dict).
    """
    if backend == "ffsim":
        sampler = ffsim.qiskit.FfsimSampler(
            default_shots=shots, norb=norb, nelec=nelec,
            global_depolarizing=depol, seed=seed,
        )
        bit_array = sampler.run([circuit]).result()[0].data.meas
        return bit_array, {"backend": backend, "noise": f"global depolarizing {depol}"}

    from qiskit import transpile
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import NoiseModel, depolarizing_error

    method = "matrix_product_state" if backend == "aer-mps" else "statevector"
    sim = AerSimulator(method=method, seed_simulator=seed)
    # optimization_level=1 and above hang on ffsim's custom LUCJ gate; level 0 lowers it
    # fine and the circuit is already hardware-shaped, so there is nothing to gain.
    transpiled = transpile(circuit, sim, optimization_level=0)

    ops = dict(transpiled.count_ops())
    one_q = [g for g in ("rz", "p", "x", "sx", "h", "u", "u1", "u2", "u3") if g in ops]
    two_q = [g for g in ("rxx", "ryy", "rzz", "cp", "cx", "cz", "swap", "xx_plus_yy")
             if g in ops]
    if depol > 0:
        noise_model = NoiseModel()
        noise_model.add_all_qubit_quantum_error(depolarizing_error(depol, 1), one_q)
        if two_q:
            noise_model.add_all_qubit_quantum_error(
                depolarizing_error(min(10 * depol, 1.0), 2), two_q
            )
        sim = AerSimulator(method=method, noise_model=noise_model, seed_simulator=seed)
    counts = sim.run(transpiled, shots=shots).result().get_counts()
    bit_array = BitArray.from_counts(counts, num_bits=2 * norb)
    return bit_array, {
        "backend": backend,
        "noise": (f"gate depolarizing 1q={depol} on {one_q}, "
                  f"2q={min(10 * depol, 1.0)} on {two_q}" if depol > 0 else "noiseless"),
        "transpiled_depth": int(transpiled.depth()),
        "transpiled_ops": {k: int(v) for k, v in ops.items()},
    }


def sampling_efficiency(bit_array, norb, nelec):
    bits = np.unpackbits(bit_array.array, axis=1, bitorder="big")[:, -2 * norb :]
    n_beta = bits[:, :norb].sum(1)
    n_alpha = bits[:, norb:].sum(1)
    return {
        "shots": int(bits.shape[0]),
        "unique_bitstrings": len({r.tobytes() for r in bit_array.array}),
        "frac_correct_particle_number": float(((n_alpha + n_beta) == sum(nelec)).mean()),
        "frac_correct_particle_number_and_sz": float(
            ((n_alpha == nelec[0]) & (n_beta == nelec[1])).mean()
        ),
    }


def sci_spin_square(sci_state, norb, nelec):
    """<S^2> of the SQD solution.

    Computed directly on the selected-CI vector via ``pyscf.fci.selected_ci``, so the
    cost tracks the *subspace* rather than the full CI space.  Embedding into the full
    vector would be simpler but allocates the whole CAS space - 5.5 GB at 44 qubits and
    77 GB at 56 - which defeats the point of a subspace method.
    """
    amps = np.asarray(sci_state.amplitudes)
    amps = amps / np.linalg.norm(amps)
    civec = fci.selected_ci._as_SCIvector(
        amps, (np.asarray(sci_state.ci_strs_a), np.asarray(sci_state.ci_strs_b))
    )
    return float(fci.selected_ci.SelectedCI().spin_square(civec, norb, nelec)[0])


def run_sqd(hcore, eri, e_nuc, bit_array, norb, nelec, spb, num_batches, max_iterations,
            seed, spin_sq=0.0, symmetrize_spin=None):
    # Spin symmetrization (the package's `open_shell=False` behaviour) unions the sampled
    # alpha and beta string sets, and is only defined when the two spin populations are
    # equal.  Default to it where it is legal, and fall back to open-shell otherwise.
    if symmetrize_spin is None:
        symmetrize_spin = nelec[0] == nelec[1]
    sci_solver = partial(solve_sci_batch, spin_sq=spin_sq, max_cycle=SCI_MAX_CYCLE)
    history = []

    def callback(results):
        best = min(results, key=lambda r: r.energy)
        history.append({
            "iteration": len(history),
            "energy": float(best.energy + e_nuc),
            "subspace_dim": int(np.prod(best.sci_state.amplitudes.shape)),
        })

    t0 = time.time()
    result = diagonalize_fermionic_hamiltonian(
        hcore, eri, bit_array,
        samples_per_batch=spb,
        norb=norb,
        nelec=nelec,
        num_batches=num_batches,
        energy_tol=1e-6,
        occupancies_tol=1e-5,
        max_iterations=max_iterations,
        sci_solver=sci_solver,
        symmetrize_spin=symmetrize_spin,
        carryover_threshold=1e-4,
        callback=callback,
        seed=seed,
    )
    occ_a, occ_b = result.orbital_occupancies
    return {
        "energy": float(result.energy + e_nuc),
        "subspace_dim": int(np.prod(result.sci_state.amplitudes.shape)),
        "occ_a": occ_a.tolist(),
        "occ_b": occ_b.tolist(),
        "spin_sq": sci_spin_square(result.sci_state, norb, nelec),
        "history": history,
        "wall_s": time.time() - t0,
    }


def oracle_curve(hcore, eri, e_nuc, norb, nelec, targets):
    """Selected-CI reference: keep the largest-|c| CASCI determinants and diagonalize in
    the product closure of their alpha/beta strings.  This is a comparison strategy, not
    an upper bound - the product closure drags in many low-weight configurations."""
    solver = fci.direct_spin1.FCI()
    solver.conv_tol = 1e-13
    solver.max_cycle = 1000
    _, civec = solver.kernel(hcore, eri, norb, nelec)
    strings_a = fci.cistring.make_strings(range(norb), nelec[0])
    strings_b = fci.cistring.make_strings(range(norb), nelec[1])
    order = np.argsort(-np.abs(civec).ravel())
    rows = []
    for n_det in targets:
        ia, ib = np.unravel_index(order[:n_det], civec.shape)
        sa, sb = np.unique(strings_a[ia]), np.unique(strings_b[ib])
        res = solve_sci((sa, sb), hcore, eri, norb, nelec, spin_sq=0.0, max_cycle=SCI_MAX_CYCLE)
        rows.append({
            "n_dets_kept": int(n_det),
            "subspace_dim": int(len(sa) * len(sb)),
            "energy": float(res.energy + e_nuc),
        })
    return rows, civec


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tag", default="sto-3g_avas0.5")
    ap.add_argument("--shots", type=int, default=100_000)
    ap.add_argument("--n-reps", type=int, default=4)
    ap.add_argument("--optimize-ansatz", action="store_true",
                    help="fit the LUCJ parameters to the t amplitudes (tutorial default)")
    ap.add_argument("--depol-sweep", type=float, nargs="+",
                    default=[0.0, 0.001, 0.003, 0.01, 0.03, 0.1])
    ap.add_argument("--depol-main", type=float, default=0.01)
    ap.add_argument("--backend", default="ffsim", choices=["ffsim", "aer", "aer-mps"],
                    help="simulator used to sample the ansatz circuit")
    ap.add_argument("--spb-sweep", type=int, nargs="+", default=[20, 50, 100, 200, 400, 800])
    ap.add_argument("--num-batches", type=int, default=3)
    ap.add_argument("--max-iterations", type=int, default=10)
    ap.add_argument("--max-full-space-dim", type=int, default=100_000,
                    help="above this CAS dimension, skip the solver-floor and oracle "
                         "diagnostics: both need the entire space, which costs hours at "
                         "19M determinants and is impossible beyond that")
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()

    s1 = load_json(RESULTS / f"stage1_{args.tag}.json")
    data = np.load(RESULTS / f"stage1_{args.tag}.npz")
    hcore, eri = data["hcore"], data["eri"]
    e_nuc = float(data["nuclear_repulsion_energy"])
    norb = s1["norb"]
    nelec = tuple(s1["nelec_active"])
    e_exact = s1["e_casci_singlet"]

    # Active spaces past exact CI have no reference (stage 1 run with --skip-casci).
    # That is the regime SQD exists for, so it must not be a crash: report energies and
    # omit the error column rather than dying on a None.
    def err(e):
        """Error string vs the exact reference, or a placeholder when there is none."""
        return f"{(e - e_exact) * 1e3:+8.4f} mHa" if e_exact is not None else "  (no ref)"
    cas_dim = s1["cas_dim"]
    ref_str = (f"{e_exact:.8f} Ha" if e_exact is not None
               else "NONE - active space is past exact CI, so SQD is unvalidated here")
    print(f"exact (CASCI) reference energy: {ref_str}   "
          f"({norb} orbitals, {2 * norb} qubits, CAS dim {cas_dim})\n", flush=True)

    circuit, ucj_op = build_circuit(
        data["t1"], data["t2"], norb, nelec, args.n_reps, args.optimize_ansatz
    )
    diag = ansatz_energy(ucj_op, hcore, eri, e_nuc, norb, nelec)
    transpiled = circuit.decompose(reps=4)
    print(f"LUCJ ansatz (n_reps={args.n_reps}): E = {diag['e_ansatz']:.8f}  "
          f"err = {err(diag['e_ansatz'])}  "
          f"PR = {diag['participation_ratio']:.2f}", flush=True)
    print(f"circuit: depth={transpiled.depth()}  "
          f"2q gates={transpiled.count_ops().get('cx', 0)}\n", flush=True)

    out = {
        "tag": args.tag,
        "config": vars(args),
        "e_casci_singlet": e_exact,
        "cas_dim": cas_dim,
        "ansatz": diag,
        "ffsim_patched": bool(ffsim_patch.PATCHED),
        "circuit_depth": int(transpiled.depth()),
        "circuit_ops": {k: int(v) for k, v in transpiled.count_ops().items()},
        "exact_occupancies": {"occ_a": s1["singlet_occ_a"], "occ_b": s1["singlet_occ_b"]},
    }

    backend_info = {}

    def sample(depol):
        bit_array, info = sample_bitstrings(
            circuit, args.shots, norb, nelec, args.backend, depol, args.seed
        )
        backend_info.update(info)
        return bit_array

    # --- experiment 1: noise sweep -------------------------------------------
    spb_fixed = args.spb_sweep[-1]
    print(f"=== noise sweep (samples_per_batch={spb_fixed}) ===", flush=True)
    rows = []
    for depol in args.depol_sweep:
        ba = sample(depol)
        eff = sampling_efficiency(ba, norb, nelec)
        res = run_sqd(hcore, eri, e_nuc, ba, norb, nelec, spb_fixed, args.num_batches,
                      args.max_iterations, args.seed)
        rows.append({"depol": depol, "efficiency": eff, **res})
        print(f"  depol={depol:<6g} valid={eff['frac_correct_particle_number_and_sz'] * 100:6.2f}%"
              f"  unique={eff['unique_bitstrings']:6d}  dim={res['subspace_dim']:7d}"
              f"  E={res['energy']:.8f}  err={err(res['energy'])}"
              f"  <S^2>={res['spin_sq']:.4f}  ({res['wall_s']:.1f}s)", flush=True)
    out["noise_sweep"] = rows

    # --- experiment 2: subspace-dimension sweep ------------------------------
    print(f"\n=== dimension sweep (depol={args.depol_main}) ===", flush=True)
    ba_main = sample(args.depol_main)
    out["main_efficiency"] = sampling_efficiency(ba_main, norb, nelec)
    rows = []
    for spb in args.spb_sweep:
        res = run_sqd(hcore, eri, e_nuc, ba_main, norb, nelec, spb, args.num_batches,
                      args.max_iterations, args.seed)
        res["samples_per_batch"] = spb
        rows.append(res)
        print(f"  spb={spb:<5d} dim={res['subspace_dim']:7d} "
              f"({res['subspace_dim'] / cas_dim * 100:5.1f}% of CAS)  "
              f"E={res['energy']:.8f}  err={err(res['energy'])}  "
              f"<S^2>={res['spin_sq']:.4f}  ({res['wall_s']:.1f}s)", flush=True)
    out["dim_sweep"] = rows

    # --- experiment 3: controls ----------------------------------------------
    print("\n=== controls ===", flush=True)
    ba_uni = generate_bit_array_uniform(args.shots, 2 * norb, rand_seed=args.seed)
    uni = run_sqd(hcore, eri, e_nuc, ba_uni, norb, nelec, spb_fixed, args.num_batches,
                  args.max_iterations, args.seed)
    uni["efficiency"] = sampling_efficiency(ba_uni, norb, nelec)
    print(f"  uniform-random: dim={uni['subspace_dim']:7d}  E={uni['energy']:.8f}  "
          f"err={err(uni['energy'])}  <S^2>={uni['spin_sq']:.4f}",
          flush=True)
    out["uniform_control"] = uni

    # Both remaining diagnostics touch the *entire* CAS space, so they are only
    # affordable on small active spaces: the solver floor diagonalizes all of it and
    # the oracle needs the exact CI vector.  At (10e,10o) that is 63,504 determinants
    # and costs seconds; at (22e,16o) it is 19,079,424 and costs hours; past that it is
    # impossible.  Guard on dimension rather than discovering it as a hung job.
    if cas_dim <= args.max_full_space_dim:
        # The selected-CI solver imposes S^2 through a Lagrange penalty and runs for a
        # bounded number of Davidson cycles, so even diagonalizing the *entire* CAS space
        # through this code path leaves a residual.  That residual is the floor any SQD
        # number in this report is measured against.
        all_strings = fci.cistring.make_strings(range(norb), nelec[0])
        floor = solve_sci((all_strings, all_strings), hcore, eri, norb, nelec,
                          spin_sq=0.0, max_cycle=SCI_MAX_CYCLE)
        out["solver_floor"] = {
            "energy": float(floor.energy + e_nuc),
            "error_mha": (None if e_exact is None
                          else float((floor.energy + e_nuc - e_exact) * 1e3)),
            "subspace_dim": int(cas_dim),
        }
        print(f"  solver floor (full CAS, same solver): E={floor.energy + e_nuc:.8f}  "
              f"err={err(floor.energy + e_nuc)}", flush=True)

        oracle, civec = oracle_curve(hcore, eri, e_nuc, norb, nelec,
                                     [25, 50, 100, 200, 400, 800])
        for row in oracle:
            print(f"  oracle top-{row['n_dets_kept']:<4d} dim={row['subspace_dim']:7d}  "
                  f"E={row['energy']:.8f}  err={err(row['energy'])}",
                  flush=True)
        out["oracle"] = oracle

        weights = np.sort(civec.ravel() ** 2)[::-1]
        cum = np.cumsum(weights)
        out["exact_wavefunction"] = {
            "participation_ratio": float(1.0 / np.sum(weights**2)),
            "max_weight": float(weights[0]),
            "dets_for_weight": {str(t): int(np.searchsorted(cum, t) + 1)
                                for t in (0.5, 0.9, 0.99, 0.999)},
        }

    else:
        msg = (f"skipped: CAS dimension {cas_dim:,} exceeds --max-full-space-dim "
               f"{args.max_full_space_dim:,}; the solver floor and the oracle both "
               f"need the full space")
        print(f"  {msg}", flush=True)
        out["solver_floor"] = {"skipped": msg}
        out["oracle"] = []
        out["exact_wavefunction"] = {"skipped": msg}

    out["backend_info"] = backend_info
    suffix = "_opt" if args.optimize_ansatz else ""
    if args.backend != "ffsim":
        suffix += f"_{args.backend}"
    save_json(RESULTS / f"stage2_{args.tag}{suffix}.json", out)
    print(f"\nwrote results/stage2_{args.tag}{suffix}.json", flush=True)


if __name__ == "__main__":
    main()

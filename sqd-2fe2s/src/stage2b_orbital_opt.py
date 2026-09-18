"""Stage 2b - improve the SQD estimate with orbital optimization.

Port of the package guide `docs/guides/optimize_orbitals.ipynb`.

Stage 2 shows that in the canonical AVAS basis SQD only reaches the exact energy once
the sampled subspace covers essentially the whole CAS space, i.e. it is exact but not
compressive. Orbital optimization is the package's answer to exactly that: the subspace
is *frozen* to the configurations SQD found and only the orbital basis is varied, so any
energy drop is compression rather than a bigger subspace.

On this system it does not help much - it recovers a few percent of the gap - which is
reported as-is in RESULTS.md.

For each subspace size we report the energy before and after optimization, so the two
numbers are directly comparable at identical subspace dimension.
"""

from __future__ import annotations

import argparse
import pathlib
import sys
import time

import ffsim
import numpy as np
from pyscf import fci

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import ffsim_patch  # noqa: E402,F401
from common import RESULTS, load_json, save_json  # noqa: E402
from stage2_sqd import SCI_MAX_CYCLE, build_circuit, sampling_efficiency  # noqa: E402
from qiskit_addon_sqd.fermion import diagonalize_fermionic_hamiltonian, solve_sci_batch  # noqa: E402
from functools import partial  # noqa: E402


def diagonalize_fixed(hamiltonian, ci_strings, norb, nelec, spin_sq):
    """Diagonalize in a fixed configuration subspace and return (energy, rdm)."""
    myci = fci.selected_ci.SelectedCI()
    myci = fci.addons.fix_spin_(myci, ss=spin_sq)
    _, amplitudes = fci.selected_ci.kernel_fixed_space(
        myci,
        hamiltonian.one_body_tensor,
        hamiltonian.two_body_tensor,
        norb,
        nelec,
        ci_strs=ci_strings,
        max_cycle=SCI_MAX_CYCLE,
    )
    dm1, dm2 = myci.make_rdm12(amplitudes, norb, nelec)
    rdm = ffsim.ReducedDensityMatrix(dm1, dm2)
    return float(rdm.expectation(hamiltonian).real), rdm


def optimize_orbitals(hcore, eri, e_nuc, ci_strings, norb, nelec, spin_sq, num_iters):
    hamiltonian = ffsim.MolecularHamiltonian(hcore, eri, constant=e_nuc)
    trace = []
    for _ in range(num_iters):
        energy, rdm = diagonalize_fixed(hamiltonian, ci_strings, norb, nelec, spin_sq)
        trace.append(energy)
        # optimize_orbitals returns U minimizing rdm.rotated(U).expectation(H),
        # equivalently rdm.expectation(H.rotated(U^dagger)), so rotate by U^dagger.
        rotation = ffsim.optimize_orbitals(rdm, hamiltonian)
        hamiltonian = hamiltonian.rotated(rotation.T.conj())
    energy, _ = diagonalize_fixed(hamiltonian, ci_strings, norb, nelec, spin_sq)
    trace.append(energy)
    return trace


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tag", default="sto-3g_avas0.5")
    ap.add_argument("--shots", type=int, default=100_000)
    ap.add_argument("--n-reps", type=int, default=4)
    ap.add_argument("--depol", type=float, default=0.01)
    ap.add_argument("--spb-sweep", type=int, nargs="+", default=[20, 50, 100, 200])
    ap.add_argument("--num-batches", type=int, default=3)
    ap.add_argument("--max-iterations", type=int, default=10)
    ap.add_argument("--oo-iters", type=int, default=10)
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()

    s1 = load_json(RESULTS / f"stage1_{args.tag}.json")
    data = np.load(RESULTS / f"stage1_{args.tag}.npz")
    hcore, eri = data["hcore"], data["eri"]
    e_nuc = float(data["nuclear_repulsion_energy"])
    norb = s1["norb"]
    nelec = tuple(s1["nelec_active"])
    e_exact = s1["e_casci_singlet"]
    cas_dim = s1["cas_dim"]
    spin_sq = 0.0

    print(f"exact (CASCI) reference energy: {e_exact:.8f} Ha  (CAS dim {cas_dim})\n",
          flush=True)

    circuit, _ = build_circuit(data["t1"], data["t2"], norb, nelec, args.n_reps, False)
    bit_array = ffsim.qiskit.FfsimSampler(
        default_shots=args.shots, norb=norb, nelec=nelec,
        global_depolarizing=args.depol, seed=args.seed,
    ).run([circuit]).result()[0].data.meas
    eff = sampling_efficiency(bit_array, norb, nelec)
    print(f"samples: {eff['shots']} shots, "
          f"{eff['frac_correct_particle_number_and_sz'] * 100:.2f}% valid, "
          f"{eff['unique_bitstrings']} unique\n", flush=True)

    sci_solver = partial(solve_sci_batch, spin_sq=spin_sq, max_cycle=SCI_MAX_CYCLE)
    rows = []
    for spb in args.spb_sweep:
        t0 = time.time()
        result = diagonalize_fermionic_hamiltonian(
            hcore, eri, bit_array,
            samples_per_batch=spb, norb=norb, nelec=nelec,
            num_batches=args.num_batches,
            energy_tol=1e-6, occupancies_tol=1e-5,
            max_iterations=args.max_iterations,
            sci_solver=sci_solver, symmetrize_spin=True,
            carryover_threshold=1e-4, seed=args.seed,
        )
        ci_strings = (result.sci_state.ci_strs_a, result.sci_state.ci_strs_b)
        dim = int(np.prod(result.sci_state.amplitudes.shape))
        e_before = float(result.energy + e_nuc)

        trace = optimize_orbitals(hcore, eri, e_nuc, ci_strings, norb, nelec,
                                  spin_sq, args.oo_iters)
        e_after = min(trace)
        rows.append({
            "samples_per_batch": spb,
            "subspace_dim": dim,
            "energy_before_oo": e_before,
            "energy_after_oo": e_after,
            "oo_trace": trace,
            "wall_s": time.time() - t0,
        })
        print(f"  spb={spb:<4d} dim={dim:7d} ({dim / cas_dim * 100:5.1f}% of CAS)  "
              f"before OO {e_before:.8f} ({(e_before - e_exact) * 1e3:+8.3f} mHa)  ->  "
              f"after OO {e_after:.8f} ({(e_after - e_exact) * 1e3:+8.3f} mHa)  "
              f"({time.time() - t0:.0f}s)", flush=True)

    save_json(RESULTS / f"stage2b_{args.tag}.json", {
        "tag": args.tag,
        "config": vars(args),
        "e_casci_singlet": e_exact,
        "cas_dim": cas_dim,
        "efficiency": eff,
        "rows": rows,
    })
    print(f"\nwrote results/stage2b_{args.tag}.json", flush=True)


if __name__ == "__main__":
    main()

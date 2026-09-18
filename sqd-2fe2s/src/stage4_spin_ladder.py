"""Stage 4 - the exchange coupling J from an SQD spin ladder.

This is the run that needs a cluster, and the one that produces a number the
iron-sulfur literature actually quotes.

Stages 1-3 validate SQD against exact CASCI on the (10e,10o) Fe 3d space.  That space
has no bridging-sulfur 3p orbitals, so it has no superexchange pathway and J comes out
roughly eight times too small.  Fixing that means a bigger active space, and past
~(26e,18o) CASCI stops being affordable - which is exactly the regime SQD is for.

So: run SQD once per Sz sector and fit the Heisenberg form to the resulting ladder.
The lowest state in the Sz = m sector is the lowest state with S >= m, so for a clean
ladder it is the S = m state; we target S^2 = m(m+1) in the solver and check what comes
back.  Fitting

    E(S) = J S(S+1) + const        (H = 2 J S1.S2, J > 0 antiferromagnetic)

gives J directly comparable to DMRG (236 cm^-1), BS-DFT (310) and experiment (148+/-16).

With --exact, the same ladder is also computed by CASCI so the two can be compared
head to head; only do that where the CAS dimension allows it.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time

import ffsim
import numpy as np
from qiskit import QuantumCircuit, QuantumRegister

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import ffsim_patch  # noqa: E402,F401
from common import HARTREE2CM, RESULTS, load_json, save_json  # noqa: E402
from stage1_classical import (  # noqa: E402
    J_REFERENCE_CM1,
    fit_heisenberg,
    lowest_in_sector,
)
from stage2_sqd import (  # noqa: E402
    SCI_MAX_CYCLE,
    run_sqd,
    sample_bitstrings,
    sampling_efficiency,
)


def build_circuit_for_sector(t1, t2, norb, nelec, n_reps, optimize):
    """LUCJ circuit with a reference in the requested Sz sector.

    The UCJ operator stays spin-balanced (it is built from the closed-shell CCSD
    amplitudes); only the Hartree-Fock reference changes, which is all that is needed
    to move between Sz sectors.  Reusing one set of amplitudes across sectors is an
    approximation, but SQD only needs the circuit as a sampler.
    """
    pairs_aa = [(p, p + 1) for p in range(norb - 1)]
    ucj_op = ffsim.UCJOpSpinBalanced.from_t_amplitudes(
        t2=t2, t1=t1, n_reps=n_reps, interaction_pairs=(pairs_aa, None), optimize=optimize
    )
    qubits = QuantumRegister(2 * norb, name="q")
    circuit = QuantumCircuit(qubits)
    circuit.append(ffsim.qiskit.PrepareHartreeFockJW(norb, nelec), qubits)
    circuit.append(ffsim.qiskit.UCJOpSpinBalancedJW(ucj_op), qubits)
    circuit.measure_all()
    return circuit


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tag", required=True,
                    help="stage-1 tag, e.g. sto-3g_fe3d or sto-3g_fe3d+brs3p")
    ap.add_argument("--shots", type=int, default=200_000)
    ap.add_argument("--n-reps", type=int, default=4)
    ap.add_argument("--optimize-ansatz", action="store_true", default=True)
    ap.add_argument("--no-optimize-ansatz", dest="optimize_ansatz", action="store_false")
    ap.add_argument("--backend", default="ffsim", choices=["ffsim", "aer", "aer-mps"])
    ap.add_argument("--depol", type=float, default=0.0)
    ap.add_argument("--samples-per-batch", type=int, default=2000)
    ap.add_argument("--num-batches", type=int, default=3)
    ap.add_argument("--max-iterations", type=int, default=10)
    ap.add_argument("--s-max", type=int, default=5)
    ap.add_argument("--spin-tol", type=float, default=0.05,
                    help="how far <S^2> may sit from S(S+1) before a sector is dropped "
                         "from the J fit; the selected-CI spin penalty is only "
                         "approximate, so this must be looser than for exact CASCI")
    ap.add_argument("--exact", action="store_true",
                    help="also compute each sector by CASCI (only where affordable)")
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()

    s1 = load_json(RESULTS / f"stage1_{args.tag}.json")
    data = np.load(RESULTS / f"stage1_{args.tag}.npz")
    hcore, eri = data["hcore"], data["eri"]
    e_nuc = float(data["nuclear_repulsion_energy"])
    norb = s1["norb"]
    n_tot = sum(s1["nelec_active"])

    print(f"active space: ({n_tot}e, {norb}o)   {2 * norb} qubits   "
          f"CAS dim {s1['cas_dim']:,}", flush=True)
    print(f"ansatz: LUCJ n_reps={args.n_reps} optimize={args.optimize_ansatz}, "
          f"backend={args.backend}, depol={args.depol}, shots={args.shots}", flush=True)
    print(f"SQD: samples_per_batch={args.samples_per_batch}, "
          f"num_batches={args.num_batches}, max_iterations={args.max_iterations}\n",
          flush=True)

    sqd_states, exact_states, rows = [], [], []
    for m in range(args.s_max + 1):
        na, nb = (n_tot + 2 * m) // 2, (n_tot - 2 * m) // 2
        if nb < 0 or na > norb:
            print(f"Sz={m}: not representable in this space, stopping", flush=True)
            break
        nelec = (na, nb)
        t0 = time.time()

        circuit = build_circuit_for_sector(
            data["t1"], data["t2"], norb, nelec, args.n_reps, args.optimize_ansatz
        )
        bit_array, backend_info = sample_bitstrings(
            circuit, args.shots, norb, nelec, args.backend, args.depol, args.seed
        )
        eff = sampling_efficiency(bit_array, norb, nelec)
        t_sample = time.time() - t0

        res = run_sqd(
            hcore, eri, e_nuc, bit_array, norb, nelec,
            args.samples_per_batch, args.num_batches, args.max_iterations,
            args.seed, spin_sq=float(m * (m + 1)),
        )
        sqd_states.append({"e_cas": res["energy"] - e_nuc, "spin_sq": res["spin_sq"],
                           "S": float(m)})
        row = {
            "sz": m, "nelec": [na, nb],
            "e_sqd": res["energy"], "spin_sq_sqd": res["spin_sq"],
            "subspace_dim": res["subspace_dim"],
            "efficiency": eff, "t_sample_s": t_sample, "t_sqd_s": res["wall_s"],
        }
        line = (f"Sz={m} ({na},{nb})  SQD E={res['energy']:.8f}  "
                f"<S^2>={res['spin_sq']:7.4f}  dim={res['subspace_dim']:,}  "
                f"sample {t_sample:.0f}s + sqd {res['wall_s']:.0f}s")

        if args.exact:
            e_ex, s_ex, s_val = lowest_in_sector(
                hcore, eri, norb, nelec, nroots=(3 if m == 0 else 1)
            )
            exact_states.append({"e_cas": e_ex, "spin_sq": s_ex, "S": s_val})
            row["e_exact"] = e_ex + e_nuc
            row["spin_sq_exact"] = s_ex
            row["error_mha"] = (res["energy"] - (e_ex + e_nuc)) * 1e3
            line += f"  | exact {e_ex + e_nuc:.8f}  err={row['error_mha']:+.3f} mHa"

        print(line, flush=True)
        rows.append(row)

    out = {
        "tag": args.tag,
        "config": vars(args),
        "norb": norb,
        "nelec_total": n_tot,
        "n_qubits": 2 * norb,
        "cas_dim": s1["cas_dim"],
        "nuclear_repulsion_energy": e_nuc,
        "rows": rows,
        "backend_info": backend_info,
        "reference_J_cm1": J_REFERENCE_CM1,
    }

    print("\n--- Heisenberg fit, H = 2 J S1.S2, E(S) = J S(S+1) ---", flush=True)
    fit_sqd = fit_heisenberg(sqd_states, spin_tol=args.spin_tol)
    out["heisenberg_sqd"] = fit_sqd
    if fit_sqd.get("ok"):
        print(f"SQD   J = {fit_sqd['J_cm1']:8.1f} cm^-1  "
              f"(residual {fit_sqd['max_residual_hartree']:.1e} Ha)", flush=True)
    else:
        print(f"SQD   fit failed: {fit_sqd['reason']}", flush=True)
    if exact_states:
        fit_ex = fit_heisenberg(exact_states)
        out["heisenberg_exact"] = fit_ex
        if fit_ex.get("ok"):
            print(f"CASCI J = {fit_ex['J_cm1']:8.1f} cm^-1", flush=True)
    for name, ref in J_REFERENCE_CM1.items():
        print(f"  reference: {name:42s} {ref:6.0f} cm^-1", flush=True)

    save_json(RESULTS / f"stage4_{args.tag}.json", out)
    print(f"\nwrote results/stage4_{args.tag}.json", flush=True)


if __name__ == "__main__":
    main()

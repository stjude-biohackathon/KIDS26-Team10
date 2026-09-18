"""Stage 5 - runtime scaling: exact CASCI vs SQD, from small to large active spaces.

This is the crossover measurement. The two methods scale completely differently:

  * **CASCI** diagonalises the *whole* space, so its cost tracks the number of
    determinants, which grows combinatorially with the number of active orbitals.
  * **SQD** diagonalises only the subspace its samples span. At fixed
    `samples_per_batch` that subspace size is roughly constant, so the cost is
    approximately flat in the active-space size - the growth shows up in circuit
    simulation, not in the eigensolver.

Plotting both against orbital count therefore shows where exact diagonalisation stops
being affordable and SQD keeps running. That crossover is the whole argument for the
method, and it is cheap to measure directly.

**Nested ladder.** We take the real (34e,22o) Fe/S Hamiltonian from stage 1 and truncate
it to the first `norb` active orbitals, at half filling. Every point is therefore a
genuine molecular Hamiltonian with realistic integral structure and sparsity, just a
smaller slice of it. The truncated *energies* are not physically meaningful - only the
largest slice is - so this is reported strictly as a **timing study**, never as
chemistry. The physically meaningful points are the real AVAS ladder in stage 1.

Timing is wall clock on one machine, single process, so the comparison between the two
methods at a given size is fair even if absolute numbers are machine-dependent.
"""

from __future__ import annotations

import argparse
import pathlib
import sys
import time

import ffsim
import numpy as np
import pyscf.fci
from scipy.special import comb

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import ffsim_patch  # noqa: E402,F401
from common import RESULTS, load_json, save_json  # noqa: E402
from stage2_sqd import run_sqd, sample_bitstrings  # noqa: E402
from stage4_spin_ladder import build_circuit_for_sector  # noqa: E402


def truncate(hcore, eri, norb):
    """Nested slice of the Hamiltonian: keep the first `norb` active orbitals."""
    return np.ascontiguousarray(hcore[:norb, :norb]), np.ascontiguousarray(
        eri[:norb, :norb, :norb, :norb]
    )


def truncate_amplitudes(t1, t2, norb, nelec):
    """Slice the CCSD amplitudes to match a truncated active space.

    t1 is (nocc, nvir) and t2 is (nocc, nocc, nvir, nvir) in the *full* space; the
    truncated space has its own occupied/virtual split, so we take the leading block of
    each and zero-pad if the source block is smaller.
    """
    nocc, nvir = nelec[0], norb - nelec[0]
    t1_new = np.zeros((nocc, nvir))
    t2_new = np.zeros((nocc, nocc, nvir, nvir))
    o = min(nocc, t1.shape[0])
    v = min(nvir, t1.shape[1])
    t1_new[:o, :v] = t1[:o, :v]
    t2_new[:o, :o, :v, :v] = t2[:o, :o, :v, :v]
    return t1_new, t2_new


def time_casci(hcore, eri, norb, nelec, max_dim, max_seconds):
    """Lowest root by exact diagonalisation, or None if it is out of budget."""
    dim = int(comb(norb, nelec[0])) * int(comb(norb, nelec[1]))
    if dim > max_dim:
        return {"skipped": f"CAS dimension {dim} exceeds --max-casci-dim {max_dim}"}
    solver = pyscf.fci.direct_spin1.FCI()
    solver.nroots = 1
    solver.conv_tol = 1e-10
    solver.max_cycle = 1000
    solver.max_space = 30
    t0 = time.time()
    energy, _ = solver.kernel(hcore, eri, norb, nelec)
    wall = time.time() - t0
    return {"energy": float(np.atleast_1d(energy)[0]), "wall_s": wall}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source-tag", default="sto-3g_fe3d+brs3p+2terms3p",
                    help="stage-1 run whose Hamiltonian is sliced (needs norb >= max)")
    ap.add_argument("--norb-list", type=int, nargs="+", default=[6, 8, 10, 12, 14])
    ap.add_argument("--shots", type=int, default=50_000)
    ap.add_argument("--n-reps", type=int, default=2)
    ap.add_argument("--samples-per-batch", type=int, default=300)
    ap.add_argument("--num-batches", type=int, default=2)
    ap.add_argument("--max-iterations", type=int, default=3)
    ap.add_argument("--no-optimize-ansatz", dest="optimize_ansatz",
                    action="store_false", default=True,
                    help="by default the ansatz parameters are fitted, which is what "
                         "makes the sampled subspace fill to a comparable size at every "
                         "orbital count - without it SQD looks artificially fast "
                         "because its subspace collapses")
    ap.add_argument("--max-casci-dim", type=int, default=20_000_000,
                    help="refuse exact diagonalisation above this many determinants")
    ap.add_argument("--max-casci-seconds", type=float, default=3600.0)
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()

    s1 = load_json(RESULTS / f"stage1_{args.source_tag}.json")
    data = np.load(RESULTS / f"stage1_{args.source_tag}.npz")
    hcore_full, eri_full = data["hcore"], data["eri"]
    norb_max = s1["norb"]
    print(f"source Hamiltonian: {args.source_tag}  ({s1['norb']} orbitals)\n", flush=True)

    if max(args.norb_list) > norb_max:
        raise SystemExit(f"--norb-list asks for {max(args.norb_list)} orbitals but the "
                         f"source has only {norb_max}")

    print(f"{'norb':>5} {'nelec':>8} {'qubits':>7} {'CAS dets':>16} "
          f"{'CASCI (s)':>11} {'SQD sample':>11} {'SQD solve':>10} "
          f"{'SQD dim':>9}", flush=True)
    print("-" * 84, flush=True)

    rows = []
    for norb in args.norb_list:
        k = norb // 2
        nelec = (k, k)
        dim = int(comb(norb, k)) ** 2
        hcore, eri = truncate(hcore_full, eri_full, norb)

        cas = time_casci(hcore, eri, norb, nelec,
                         args.max_casci_dim, args.max_casci_seconds)

        t1, t2 = truncate_amplitudes(data["t1"], data["t2"], norb, nelec)
        t0 = time.time()
        circuit = build_circuit_for_sector(
            t1, t2, norb, nelec, args.n_reps, args.optimize_ansatz
        )
        bit_array, _ = sample_bitstrings(
            circuit, args.shots, norb, nelec, "ffsim", 0.0, args.seed
        )
        t_sample = time.time() - t0

        sqd = run_sqd(hcore, eri, 0.0, bit_array, norb, nelec,
                      args.samples_per_batch, args.num_batches,
                      args.max_iterations, args.seed, spin_sq=0.0)

        cas_str = f"{cas['wall_s']:11.2f}" if "wall_s" in cas else f"{'skipped':>11}"
        print(f"{norb:5d} {str(nelec):>8} {2 * norb:7d} {dim:16,} {cas_str} "
              f"{t_sample:11.2f} {sqd['wall_s']:10.2f} {sqd['subspace_dim']:9,}",
              flush=True)

        rows.append({
            "norb": norb, "nelec": list(nelec), "n_qubits": 2 * norb,
            "cas_dim": dim,
            "casci": cas,
            "sqd_sample_s": t_sample,
            "sqd_solve_s": sqd["wall_s"],
            "sqd_subspace_dim": sqd["subspace_dim"],
            "sqd_energy": sqd["energy"],
            "circuit_depth": int(circuit.decompose(reps=4).depth()),
        })

    save_json(RESULTS / "stage5_scaling.json", {
        "source_tag": args.source_tag,
        "config": vars(args),
        "note": "TIMING STUDY ONLY - truncated active spaces are not physically "
                "meaningful; energies are reported for reproducibility, not chemistry",
        "rows": rows,
    })
    print("\nwrote results/stage5_scaling.json", flush=True)


if __name__ == "__main__":
    main()

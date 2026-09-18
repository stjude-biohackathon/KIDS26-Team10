"""Consolidate every real-QPU run into one file: results/hardware_all.json.

Run this after any stage7 job (or replay) finishes:

    python3 src/collect_hardware.py

Why this exists: the 20-qubit run predates stage7's --tag/--suffix flags, so its result
file (`stage7_hardware_result.json`) was later overwritten by the 32-qubit run.  The
20-qubit numbers are therefore pinned here from its run log, which is kept in
`results/qpu_run_20q.log`; its raw samples are still on disk if anyone wants to
re-derive them.  Every other run is read straight from its stage7 JSON.
"""
from __future__ import annotations

import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from common import RESULTS, save_json  # noqa: E402

# qubits, active space, plain/mitigated, stage7 file (None = pinned from log)
RUNS = [
    (20, "(10e,10o)", False, None),
    (32, "(22e,16o)", False, "stage7_hardware_result.json"),
    (32, "(22e,16o)", True, "stage7_hardware_sto-3g_fe3d+brs3p_mitigated.json"),
    (36, "(26e,18o)", False, "stage7_hardware_sto-3g_fe3d+brs3p+fe4s.json"),
    (36, "(26e,18o)", True, "stage7_hardware_sto-3g_fe3d+brs3p+fe4s_mitigated.json"),
]

# The 20-qubit run, from results/qpu_run_20q.log (job dameitg2fm4c73f2t3f0).
PINNED_20Q = {
    "backend": "ibm_fez", "job_id": "dameitg2fm4c73f2t3f0", "shots": 100000,
    "transpiled_depth": 1677, "transpiled_2q": 2031,
    "usage": "{'quantum_seconds': None}",
    "e_casci_exact": -5013.64906670,
    "raw_efficiency": {"shots": 100000, "unique_bitstrings": 95300,
                       "frac_correct_particle_number": 0.1722,
                       "frac_correct_particle_number_and_sz": 0.0590},
    "sqd_energy": -5013.64906670, "sqd_error_mha": 0.0,
    "sqd_subspace_dim": 63504, "sqd_spin_sq": 0.0,
}

# total determinants in each active space, for the coverage column
SPACE_DIM = {20: 63504, 32: 19079424, 36: 73410624}
# extrapolated singlet at (26e,18o) - NOT exact; from the four exact spin sectors
EXTRAP_36Q = -5013.67603364


def usage_seconds(s):
    try:
        v = json.loads(str(s).replace("'", '"').replace("None", "null"))
        return v.get("quantum_seconds")
    except Exception:
        return None


def main():
    runs, missing = [], []
    for qb, space, mit, fname in RUNS:
        if fname is None:
            d = PINNED_20Q
        else:
            p = RESULTS / fname
            if not p.exists():
                missing.append(fname)
                continue
            d = json.load(open(p))
        eff = d["raw_efficiency"]
        ref = d.get("e_casci_exact")
        err = d.get("sqd_error_mha")
        if ref is None and d.get("sqd_energy") is not None:
            ref, err = EXTRAP_36Q, (d["sqd_energy"] - EXTRAP_36Q) * 1e3
        dim = d.get("sqd_subspace_dim")
        runs.append({
            "qubits": qb, "active_space": space, "job_id": d["job_id"],
            "mitigation": ("dynamical decoupling XpXm + gate/measure twirling"
                           if mit else None),
            "shots": d["shots"],
            "transpiled_depth": d["transpiled_depth"],
            "two_qubit_gates": d["transpiled_2q"],
            "measured_quantum_seconds": usage_seconds(d.get("usage")),
            "frac_valid_N": eff["frac_correct_particle_number"],
            "frac_valid_N_and_Sz": eff["frac_correct_particle_number_and_sz"],
            "unique_bitstrings": eff["unique_bitstrings"],
            "energy": d.get("sqd_energy"), "reference": ref,
            "error_mha": None if err is None or err != err else err,
            "spin_sq": d.get("sqd_spin_sq"),
            "subspace_dim": dim,
            "pct_of_space": None if dim is None else 100 * dim / SPACE_DIM[qb],
            "reference_kind": ("exact CASCI" if d.get("e_casci_exact") is not None
                               else "EXTRAPOLATED from four exact spin sectors"),
        })
    secs = [r["measured_quantum_seconds"] for r in runs
            if r["measured_quantum_seconds"]]
    out = {
        "backend": "ibm_fez (156-qubit Heron)",
        "shots_per_job": 100000,
        "n_jobs_completed": len(runs),
        "measured_quantum_seconds_reported": secs,
        "measured_quantum_seconds_total_reported": sum(secs),
        "usage_note": (
            "IBM reported quantum_seconds for three of the jobs; the mean is "
            f"{sum(secs)/len(secs):.1f} s per 100,000 shots, i.e. ~"
            f"{1e6*sum(secs)/len(secs)/100000:.0f} us per shot. Our pre-run estimates "
            "used Qiskit's estimate_duration, which counts gate time only and omits "
            "inter-shot reset/thermalisation, so they were ~6x too low. Estimated "
            "total spend across all five completed jobs is ~250 s of a 600 s budget."),
        "cancelled": {
            "job_id": "damf1k02fm4c73f2trag", "shots": 1000000, "qubits": 32,
            "reason": ("IBM's own queue estimate was ~13 min, over the 10-minute "
                       "allocation; cancelled before execution")},
        "runs": runs,
    }
    save_json(RESULTS / "hardware_all.json", out)
    for r in runs:
        e = "     n/a" if r["error_mha"] is None else f"{r['error_mha']:+9.4f}"
        print(f"  {r['qubits']:>2}q {'MIT ' if r['mitigation'] else'plain'} "
              f"{r['two_qubit_gates']:>5} 2q  N={r['frac_valid_N']*100:5.2f}%  "
              f"N+Sz={r['frac_valid_N_and_Sz']*100:5.2f}%  err={e} mHa  {r['job_id']}")
    if missing:
        print(f"  NOT YET ON DISK: {missing}")
    print(f"wrote {RESULTS / 'hardware_all.json'}")


if __name__ == "__main__":
    main()

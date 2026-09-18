"""SQD on REAL IBM Quantum hardware counts, compared against the exact
CASCI reference, the simulator-based dense-LUCJ SQD run, and the matched-
dimension random baseline (all three already computed in
src.sqd_benchmark).

Pipeline
--------
1. Load the raw measurement `BitArray` from the real `ibm_kingston` job
   (src.hardware_run: job_id saved in data/processed/hardware_job.json,
   counts fetched losslessly as packed bits in data/processed/
   hardware_counts.npz -- NOT reconstructed from the lossy hex-string
   `get_counts()` dict).
2. Run the IDENTICAL SQD self-consistent configuration-recovery loop
   (`fermion.diagonalize_fermionic_hamiltonian`, same samples_per_batch/
   num_batches/max_iterations as src.sqd_benchmark) on those hardware
   counts against the same frozen restricted (10e,10o) Hamiltonian.
3. Report E_SQD(hardware) alongside E_CASCI (exact), E_SQD(dense-LUCJ,
   simulator) and E_baseline (matched-shot-budget uniform random), plus
   raw-shot sampling diagnostics (unique bitstrings, Shannon entropy) for
   all three samplers side by side.

This is a genuinely different circuit from the simulator LUCJ benchmark:
src.lucj_ansatz_hw's hardware-native ("Local" UCJ, linear interaction_pairs)
ansatz, optimized separately (537.5 mHa gap to CASCI vs the dense ansatz's
own optimum) specifically because it is the one that transpiles shallow
enough (338 CZ gates / depth 235 on ibm_kingston, vs 1658 CZ / depth 1188
for the dense ansatz) to have a physically meaningful (~50% at median gate
error, per src.hardware_run's docstring) expected fidelity on today's
hardware.
"""

import json

import numpy as np

from src.restricted_hamiltonian import load_restricted_hamiltonian
from src.sqd_benchmark import run_sqd, sampling_diagnostics, SQD_BENCHMARK_PATH
from src.hardware_run import load_hardware_bit_array, JOB_META_PATH

HARDWARE_SQD_PATH = "data/processed/hardware_sqd_benchmark.npz"


if __name__ == "__main__":
    d = load_restricted_hamiltonian()
    h1e, eri, ecore = d["h1e"], d["eri"], d["ecore"]
    ncas, nelecas = d["ncas"], d["nelecas"]
    e_casci = d["e_casci"]

    with open(JOB_META_PATH) as f:
        job_meta = json.load(f)
    print(f"real-hardware job: {job_meta['job_id']} on {job_meta['backend']}, "
          f"{job_meta['shots']} shots, {job_meta['n_2q_gates']} 2-qubit gates, "
          f"depth {job_meta['depth']}")

    bit_array_hw = load_hardware_bit_array()
    n_shots_hw = bit_array_hw.num_shots
    diag_hw = sampling_diagnostics(bit_array_hw, n_shots_hw)
    print(f"raw shot diversity (hardware): {diag_hw['n_unique']} unique bitstrings / "
          f"{n_shots_hw} shots, Shannon entropy = {diag_hw['entropy_bits']:.2f} bits "
          f"(max possible = {np.log2(n_shots_hw):.2f} bits)")

    print(f"\nrunning SQD self-consistent configuration recovery on "
          f"REAL HARDWARE samples ...")
    sci_hw, dt_hw, dim_hw = run_sqd(h1e, eri, bit_array_hw, ncas, nelecas)
    e_sqd_hw = sci_hw.energy + ecore
    print(f"  E_SQD(hardware) = {e_sqd_hw:.8f} Ha  (subspace dim {dim_hw})  "
          f"({dt_hw:.1f} s)")
    print(f"  error vs CASCI: {(e_sqd_hw-e_casci)*1000:.4f} mHa")

    prior = np.load(SQD_BENCHMARK_PATH)
    e_sqd_sim = float(prior["e_sqd"])
    e_baseline = float(prior["e_baseline"])
    dim_sim = tuple(int(x) for x in np.atleast_1d(prior["dim_lucj"]))
    dim_rand = tuple(int(x) for x in np.atleast_1d(prior["dim_rand"]))
    n_unique_sim = int(prior["n_unique_lucj"])
    entropy_sim = float(prior["entropy_lucj"])
    n_unique_rand = int(prior["n_unique_rand"])
    entropy_rand = float(prior["entropy_rand"])

    print(f"\n{'':28}{'E (Ha)':>18}{'error (mHa)':>14}{'subspace dim':>14}"
          f"{'unique bits':>13}{'entropy (bits)':>16}")
    print(f"{'CASCI (exact)':28}{e_casci:>18.8f}{0.0:>14.4f}{'-':>14}{'-':>13}{'-':>16}")
    print(f"{'SQD (dense LUCJ, sim)':28}{e_sqd_sim:>18.8f}"
          f"{(e_sqd_sim-e_casci)*1000:>14.4f}{str(dim_sim):>14}{n_unique_sim:>13}"
          f"{entropy_sim:>16.2f}")
    print(f"{'SQD (hw-native LUCJ, HW)':28}{e_sqd_hw:>18.8f}"
          f"{(e_sqd_hw-e_casci)*1000:>14.4f}{str(dim_hw):>14}{diag_hw['n_unique']:>13}"
          f"{diag_hw['entropy_bits']:>16.2f}")
    print(f"{'random baseline':28}{e_baseline:>18.8f}"
          f"{(e_baseline-e_casci)*1000:>14.4f}{str(dim_rand):>14}{n_unique_rand:>13}"
          f"{entropy_rand:>16.2f}")

    np.savez(
        HARDWARE_SQD_PATH,
        e_casci=e_casci, e_sqd_hw=e_sqd_hw, e_sqd_sim=e_sqd_sim, e_baseline=e_baseline,
        dim_hw=np.array(dim_hw), dt_hw=dt_hw,
        n_shots_hw=n_shots_hw, n_unique_hw=diag_hw["n_unique"],
        entropy_hw=diag_hw["entropy_bits"],
        job_id=job_meta["job_id"], backend=job_meta["backend"],
        n_2q_gates=job_meta["n_2q_gates"], depth=job_meta["depth"],
        e_hf_hw=job_meta["e_hf"], e_opt_hw=job_meta["e_opt"],
    )
    print(f"\nsaved -> {HARDWARE_SQD_PATH}")

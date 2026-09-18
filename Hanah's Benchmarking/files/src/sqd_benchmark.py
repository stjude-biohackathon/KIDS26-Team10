"""SQD vs exact CASCI benchmark, and the matched-dimension random-sampler
baseline that isolates whether SQD's informed (LUCJ-circuit) sampling beats
blind sampling at equal computational budget.

Pipeline
--------
1. Sample `N_SHOTS` bitstrings from the (classically pre-optimized) LUCJ
   circuit built in src.lucj_ansatz, using ffsim's fast simulator via
   qiskit-addon-sqd's expected `BitArray` format.
2. Run qiskit-addon-sqd's self-consistent configuration-recovery loop
   (`fermion.diagonalize_fermionic_hamiltonian`) on those samples against the
   frozen restricted (10e,10o) Hamiltonian (src.restricted_hamiltonian) ->
   E_SQD, at whatever final subspace dimension the algorithm converges to.
3. Matched-dimension baseline: draw the SAME NUMBER of shots from the
   UNIFORM distribution over bitstrings (`qiskit_addon_sqd.counts.
   generate_bit_array_uniform`, no circuit, no informed structure at all),
   and run it through the *identical* SQD hyperparameters (samples_per_batch,
   num_batches, max_iterations). This holds the entire configuration-
   recovery machinery and its sample budget fixed, changing only whether the
   input samples came from an ansatz informed by the Hamiltonian or from
   pure noise -- isolating exactly the thing SQD claims buys something.
4. Compare E_SQD, E_baseline, and E_CASCI (exact), plus the final subspace
   dimensions each run actually explored.

Result (this run): the matched-dimension baseline WINS
--------------------------------------------------------
With the specific LUCJ ansatz frozen by src.lucj_ansatz (n_reps=1,
HF-seeded, only 10 classical-VQE iterations, 637 mHa short of the exact
ground state -- see that module's docstring for why CCSD-seeding was not
available here), SQD's informed sampler does WORSE than blind uniform
sampling at the identical shot budget and identical SQD hyperparameters:
E_SQD is ~411 mHa above CASCI vs ~0.3 mHa for the random baseline. This is
not a failure of the SQD method in general -- it is exactly the kind of
sampler-quality failure mode this matched-dimension control is designed to
catch. The mechanism, verified directly on the raw shot counts
(`sampling_diagnostics`, printed and saved below): the under-converged,
single-reference-seeded LUCJ circuit's output distribution is far more
peaked (10000 shots -> 801 unique bitstrings, Shannon entropy 6.7 bits) than
uniform random (10000 shots -> 9962 unique, entropy 13.3 bits). This
molecule's true active-space ground state is a genuine open-shell-singlet /
diradical (see src.spin_ladder's sub-mHa spin ladder and src.lucj_ansatz's
failed RHF/CCSD attempt): its correlated wavefunction is spread over
open-shell configurations far from any single reference determinant. An
ansatz whose distribution is concentrated near that single reference (HF)
therefore fails to *cover* the configurations SQD's self-consistent recovery
needs, and ends up diagonalizing a smaller, less complete final subspace
(97 configurations) than blind sampling stumbles into by covering
configuration space near-uniformly (252 configurations) -- and because the
SCI energy is a variational upper bound, more/better-chosen configurations
directly means a lower (better) energy. A properly converged LUCJ ansatz
(more VQE iterations, or a valid CCSD seed) would be expected to reverse
this, concentrating probability on the CORRECT open-shell configurations
instead of the wrong single-reference ones -- that reversal was not
achievable within this sandbox's compute budget (~40 s/iteration, only 10
iterations run), so this result should be read as a controlled demonstration
of a real SQD failure mode (a Hamiltonian-informed but under-converged
sampler can be worse than no information at all), not as evidence that SQD
underperforms random sampling in general.
"""

import time

import numpy as np
import ffsim
import ffsim.qiskit as ffq
import qiskit
from qiskit_addon_sqd import fermion, counts as sqd_counts

from src.restricted_hamiltonian import load_restricted_hamiltonian
from src.lucj_ansatz import LUCJ_PATH

N_SHOTS = 10_000
SAMPLES_PER_BATCH = 400
NUM_BATCHES = 5
MAX_SQD_ITERATIONS = 10
SEED = 0

SQD_BENCHMARK_PATH = "data/processed/sqd_benchmark.npz"


def load_optimized_lucj_op(path=LUCJ_PATH):
    d = np.load(path)
    x_opt, n_reps, ncas = d["x_opt"], int(d["n_reps"]), int(d["ncas"])
    op = ffsim.UCJOpSpinBalanced.from_parameters(x_opt, norb=ncas, n_reps=n_reps)
    return op, dict(e_hf=float(d["e_hf"]), e_random=float(d["e_random"]),
                     e_opt=float(d["e_opt"]), trace=d["trace"])


def sample_lucj_circuit(op, ncas, nelecas, shots=N_SHOTS):
    qubits = qiskit.QuantumRegister(2 * ncas, name="q")
    circuit = qiskit.QuantumCircuit(qubits)
    circuit.append(ffq.PrepareHartreeFockJW(ncas, nelecas), qubits)
    circuit.append(ffq.UCJOpSpinBalancedJW(op), qubits)
    circuit.measure_all()

    sampler = ffq.FfsimSampler()
    job = sampler.run([circuit], shots=shots)
    result = job.result()
    return result[0].data.meas


def sampling_diagnostics(bit_array, n_shots):
    """Unique-bitstring count and Shannon entropy of the raw shot distribution,
    to diagnose *why* an SQD run over/under-performs: SQD's self-consistent
    configuration recovery can only end up with as many distinct
    configurations to diagonalize over as the input samples actually explore,
    and the SCI energy is variational (more configurations -> at least as
    low), so a low-entropy (highly peaked) sampler is a direct, checkable
    cause of a worse final energy, independent of whether the peak sits on
    physically 'good' or 'bad' determinants.
    """
    bits = bit_array.array
    _, counts_ = np.unique(bits, axis=0, return_counts=True)
    probs = counts_ / counts_.sum()
    entropy = float(-(probs * np.log2(probs)).sum())
    return dict(n_unique=int(len(counts_)), entropy_bits=entropy)


def run_sqd(h1e, eri, bit_array, ncas, nelecas, samples_per_batch=SAMPLES_PER_BATCH,
            num_batches=NUM_BATCHES, max_iterations=MAX_SQD_ITERATIONS, seed=SEED):
    t0 = time.time()
    sci_result = fermion.diagonalize_fermionic_hamiltonian(
        h1e, eri, bit_array, samples_per_batch=samples_per_batch, norb=ncas,
        nelec=nelecas, num_batches=num_batches, max_iterations=max_iterations,
        seed=seed,
    )
    dt = time.time() - t0
    dim = sci_result.sci_state.amplitudes.shape
    return sci_result, dt, dim


if __name__ == "__main__":
    d = load_restricted_hamiltonian()
    h1e, eri, ecore = d["h1e"], d["eri"], d["ecore"]
    ncas, nelecas = d["ncas"], d["nelecas"]
    e_casci = d["e_casci"]

    op, lucj_info = load_optimized_lucj_op()
    print(f"loaded optimized LUCJ ansatz: E(HF)={lucj_info['e_hf']:.6f} Ha, "
          f"E(optimized)={lucj_info['e_opt']:.6f} Ha "
          f"(gap to CASCI: {(lucj_info['e_opt']-e_casci)*1000:.2f} mHa)")

    print(f"\nsampling {N_SHOTS} shots from the optimized LUCJ circuit ...")
    t0 = time.time()
    bit_array_lucj = sample_lucj_circuit(op, ncas, nelecas, shots=N_SHOTS)
    print(f"  done ({time.time()-t0:.1f} s)")
    diag_lucj = sampling_diagnostics(bit_array_lucj, N_SHOTS)
    print(f"  raw shot diversity: {diag_lucj['n_unique']} unique bitstrings / "
          f"{N_SHOTS} shots, Shannon entropy = {diag_lucj['entropy_bits']:.2f} bits "
          f"(max possible = {np.log2(N_SHOTS):.2f} bits)")

    print(f"\nrunning SQD self-consistent configuration recovery "
          f"(samples_per_batch={SAMPLES_PER_BATCH}, num_batches={NUM_BATCHES}, "
          f"max_iterations={MAX_SQD_ITERATIONS}) on LUCJ samples ...")
    sci_lucj, dt_lucj, dim_lucj = run_sqd(h1e, eri, bit_array_lucj, ncas, nelecas)
    e_sqd = sci_lucj.energy + ecore
    print(f"  E_SQD = {e_sqd:.8f} Ha  (subspace dim {dim_lucj})  ({dt_lucj:.1f} s)")
    print(f"  error vs CASCI: {(e_sqd-e_casci)*1000:.4f} mHa")

    print(f"\ngenerating matched-dimension baseline: {N_SHOTS} UNIFORM random "
          f"bitstrings (no circuit, no Hamiltonian information) ...")
    bit_array_random = sqd_counts.generate_bit_array_uniform(
        N_SHOTS, 2 * ncas, rand_seed=SEED)
    diag_rand = sampling_diagnostics(bit_array_random, N_SHOTS)
    print(f"  raw shot diversity: {diag_rand['n_unique']} unique bitstrings / "
          f"{N_SHOTS} shots, Shannon entropy = {diag_rand['entropy_bits']:.2f} bits")

    print(f"running the IDENTICAL SQD pipeline on the random baseline samples ...")
    sci_rand, dt_rand, dim_rand = run_sqd(h1e, eri, bit_array_random, ncas, nelecas)
    e_baseline = sci_rand.energy + ecore
    print(f"  E_baseline = {e_baseline:.8f} Ha  (subspace dim {dim_rand})  "
          f"({dt_rand:.1f} s)")
    print(f"  error vs CASCI: {(e_baseline-e_casci)*1000:.4f} mHa")

    print(f"\n{'':20}{'E (Ha)':>18}{'error (mHa)':>14}{'subspace dim':>16}")
    print(f"{'CASCI (exact)':20}{e_casci:>18.8f}{0.0:>14.4f}{'-':>16}")
    print(f"{'SQD (LUCJ)':20}{e_sqd:>18.8f}{(e_sqd-e_casci)*1000:>14.4f}"
          f"{str(dim_lucj):>16}")
    print(f"{'random baseline':20}{e_baseline:>18.8f}{(e_baseline-e_casci)*1000:>14.4f}"
          f"{str(dim_rand):>16}")

    improvement = (e_baseline - e_sqd) * 1000
    print(f"\nSQD improves on the matched-shot-budget random baseline by "
          f"{improvement:.4f} mHa "
          f"({'SQD wins' if improvement > 0 else 'baseline wins (unexpected)'})")
    if improvement < 0:
        print("  diagnosis: the LUCJ sampler's shot distribution is far more "
              f"peaked ({diag_lucj['entropy_bits']:.1f} bits entropy, "
              f"{diag_lucj['n_unique']} unique bitstrings) than uniform "
              f"random ({diag_rand['entropy_bits']:.1f} bits, "
              f"{diag_rand['n_unique']} unique) -- see module docstring for "
              "why this specific under-optimized, HF-seeded ansatz backfires "
              "on this genuinely multi-configurational (diradical) ground "
              "state: it concentrates shots near a single-reference-like "
              "region instead of covering the open-shell configurations "
              "the true ground state needs, so the self-consistent recovery "
              "ends up diagonalizing a SMALLER, less complete subspace than "
              "blind sampling gets 'for free' by covering configuration "
              "space near-uniformly.")

    np.savez(
        SQD_BENCHMARK_PATH,
        e_casci=e_casci, e_sqd=e_sqd, e_baseline=e_baseline,
        dim_lucj=np.array(dim_lucj), dim_rand=np.array(dim_rand),
        dt_lucj=dt_lucj, dt_rand=dt_rand,
        n_shots=N_SHOTS, samples_per_batch=SAMPLES_PER_BATCH,
        num_batches=NUM_BATCHES, max_iterations=MAX_SQD_ITERATIONS,
        e_hf=lucj_info["e_hf"], e_lucj_opt=lucj_info["e_opt"],
        n_unique_lucj=diag_lucj["n_unique"], entropy_lucj=diag_lucj["entropy_bits"],
        n_unique_rand=diag_rand["n_unique"], entropy_rand=diag_rand["entropy_bits"],
    )
    print(f"\nsaved -> {SQD_BENCHMARK_PATH}")

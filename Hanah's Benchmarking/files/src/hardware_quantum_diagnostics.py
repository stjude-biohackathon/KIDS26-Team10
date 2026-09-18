"""Genuine per-qubit / pairwise-correlation quantum diagnostics comparing the
REAL `ibm_kingston` hardware samples against the classically-simulated dense
LUCJ ansatz and a matched-budget uniform-random control -- built because a
single most-sampled bitstring (the notebook's original occupation plot) is
not informative: with only 81/10000 shots on the top bitstring, one sample
says nothing about how the full distribution behaves. Everything here is
computed from the FULL shot distribution instead.

Three diagnostics, each computed identically for all three samplers
(hardware, simulator dense-LUCJ, uniform random):

1. Shot-averaged per-qubit occupation ``<n_q>``, alpha and beta channels
   separately -- the real, statistically meaningful analogue of "the
   occupation pattern", replacing the single-bitstring plot.
2. Pairwise qubit-qubit connected correlation matrix
   ``C_ij = <n_i n_j> - <n_i><n_j>`` over all 20 qubits. The LUCJ ansatz's
   `interaction_pairs` (the specific pairs of orbitals its diagonal-Coulomb
   gates entangle) is a concrete, checkable prediction for WHERE structured
   correlation should appear. This directly tests whether that structure
   is still visible in the real-hardware statistics, or whether hardware
   noise has washed it out toward the random control's correlation matrix.
3. Excitation-order histogram: Hamming distance of each sampled bitstring
   from the Hartree-Fock reference, weighted by shot count. A sampler whose
   output still carries information from a physically-informed ansatz
   should concentrate at LOW excitation order relative to unstructured
   uniform sampling over the same 20 bits.

Reproducibility note on the simulator distribution
----------------------------------------------------
`src.sqd_benchmark.sample_lucj_circuit` does not seed `ffsim.qiskit.
FfsimSampler`, so re-sampling it here (needed because the original run's raw
bitstrings/shots were never persisted, only summary statistics) reproduces
the same PHYSICS but not a bit-identical shot record: this run found 835
unique bitstrings / 6.71 bits entropy vs. the original 878 / 6.79 bits --
close enough to draw the same qualitative conclusions, but reported here
plainly as a fresh, seeded resample rather than a byte-for-byte replay.
"""

import numpy as np
import ffsim.qiskit as ffq

from src.restricted_hamiltonian import load_restricted_hamiltonian
from src.lucj_ansatz import LUCJ_PATH
from src.sqd_benchmark import load_optimized_lucj_op, N_SHOTS, SEED
from src.hardware_run import load_hw_ansatz, build_logical_circuit, load_hardware_bit_array
from qiskit_addon_sqd import counts as sqd_counts

OUT_PATH = "data/processed/hardware_quantum_diagnostics.npz"


def qubit_ordered_bool_array(bit_array):
    """(shots, num_bits) boolean array with column q = qubit q (column 0 =
    qubit 0), matching the `reversed(bitstring)` convention used elsewhere
    in this project -- verified directly: `to_bool_array()[:, ::-1]` equals
    `[int(c) for c in reversed(get_bitstrings()[shot])]` for shot 0 of the
    real-hardware BitArray.
    """
    return bit_array.to_bool_array()[:, ::-1]


def per_qubit_occupation(bits):
    return bits.mean(axis=0)


def correlation_matrix(bits):
    bits_f = bits.astype(float)
    n = bits_f.mean(axis=0)
    cov = (bits_f.T @ bits_f) / bits_f.shape[0] - np.outer(n, n)
    return cov


def excitation_hamming_histogram(bits, hf_bits, nbits):
    hamming = (bits != hf_bits[None, :]).sum(axis=1)
    hist = np.bincount(hamming, minlength=nbits + 1)
    return hist / hist.sum()


def sample_dense_lucj(op, ncas, nelecas, shots=N_SHOTS, seed=SEED):
    """Same construction as src.sqd_benchmark.sample_lucj_circuit, but with
    an explicit seed on FfsimSampler for a reproducible diagnostic resample.
    """
    import qiskit
    qubits = qiskit.QuantumRegister(2 * ncas, name="q")
    circuit = qiskit.QuantumCircuit(qubits)
    circuit.append(ffq.PrepareHartreeFockJW(ncas, nelecas), qubits)
    circuit.append(ffq.UCJOpSpinBalancedJW(op), qubits)
    circuit.measure_all()
    sampler = ffq.FfsimSampler(seed=seed)
    job = sampler.run([circuit], shots=shots)
    return job.result()[0].data.meas


if __name__ == "__main__":
    d_ham = load_restricted_hamiltonian()
    ncas, nelecas = d_ham["ncas"], d_ham["nelecas"]
    nbits = 2 * ncas

    hf_bits = np.array([1] * nelecas[0] + [0] * (ncas - nelecas[0])
                        + [1] * nelecas[1] + [0] * (ncas - nelecas[1]), dtype=bool)

    # --- real hardware samples (already fetched, lossless) -----------------
    ba_hw = load_hardware_bit_array()
    bits_hw = qubit_ordered_bool_array(ba_hw)
    print(f"hardware: {bits_hw.shape[0]} shots loaded")

    # --- simulator dense-LUCJ resample (seeded, for this diagnostic) ------
    op_dense, _ = load_optimized_lucj_op()
    ba_sim = sample_dense_lucj(op_dense, ncas, nelecas)
    bits_sim = qubit_ordered_bool_array(ba_sim)
    print(f"simulator dense-LUCJ resample: {bits_sim.shape[0]} shots, "
          f"{len(np.unique(bits_sim, axis=0))} unique bitstrings")

    # --- matched-budget uniform random control ------------------------------
    ba_rand = sqd_counts.generate_bit_array_uniform(N_SHOTS, nbits, rand_seed=SEED)
    bits_rand = qubit_ordered_bool_array(ba_rand)
    print(f"uniform random control: {bits_rand.shape[0]} shots")

    results = {}
    for name, bits in [("hw", bits_hw), ("sim", bits_sim), ("rand", bits_rand)]:
        occ = per_qubit_occupation(bits)
        results[f"occ_a_{name}"] = occ[:ncas]
        results[f"occ_b_{name}"] = occ[ncas:]
        results[f"corr_{name}"] = correlation_matrix(bits)
        results[f"hamming_hist_{name}"] = excitation_hamming_histogram(bits, hf_bits, nbits)

    occ_a_hf = hf_bits[:ncas].astype(float)
    occ_b_hf = hf_bits[ncas:].astype(float)

    # hardware-native ansatz's own entangling structure, to overlay on the
    # hardware/simulator correlation matrices as "where structure is predicted"
    _, (pairs_aa, pairs_ab), _ = load_hw_ansatz()

    print("\nper-qubit occupation, alpha channel (orbital: hw / sim / rand / HF):")
    for p in range(ncas):
        print(f"  {p}: {results['occ_a_hw'][p]:.3f} / {results['occ_a_sim'][p]:.3f} / "
              f"{results['occ_a_rand'][p]:.3f} / {occ_a_hf[p]:.0f}")

    print("\nmean |off-diagonal correlation| (structure strength):")
    for name in ("hw", "sim", "rand"):
        c = results[f"corr_{name}"]
        off = c[~np.eye(nbits, dtype=bool)]
        print(f"  {name}: {np.abs(off).mean():.5f}")

    print("\nHamming-distance-from-HF distribution, mean +/- std:")
    for name in ("hw", "sim", "rand"):
        h = results[f"hamming_hist_{name}"]
        d = np.arange(nbits + 1)
        mean = (h * d).sum()
        std = np.sqrt((h * (d - mean) ** 2).sum())
        print(f"  {name}: mean={mean:.2f}, std={std:.2f}")

    np.savez(
        OUT_PATH,
        occ_a_hf=occ_a_hf, occ_b_hf=occ_b_hf,
        pairs_aa=np.array(pairs_aa), pairs_ab=np.array(pairs_ab),
        ncas=ncas, nelecas=np.array(nelecas), nbits=nbits,
        **results,
    )
    print(f"\nsaved -> {OUT_PATH}")

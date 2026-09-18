"""Submit the hardware-native LUCJ circuit (src.lucj_ansatz_hw) to real IBM
Quantum hardware, and later recover the sampled counts for SQD.

Credentials
-----------
This module NEVER reads or stores an API key/CRN itself. It only calls
`QiskitRuntimeService()` with no arguments, which reads whatever account
was already saved (out-of-workspace, at `~/.qiskit/qiskit-ibm-runtime.json`)
by a one-time `QiskitRuntimeService.save_account(...)` call made directly
in an interactive shell, never written to a file under `data/` or `src/`.

Backend and circuit
--------------------
Backend: `ibm_kingston` (156-qubit Heron r2, verified newest-online and
fastest-queue of the 3 backends visible to this account at submission
time). Circuit: the hardware-native ("Local" UCJ) ansatz from
src.lucj_ansatz_hw, transpiled with `ffsim.qiskit.generate_lucj_pass_manager`
(heavy-hex-aware layout, not generic `generate_preset_pass_manager`) --
verified directly (see that module's docstring) to produce 338 two-qubit
(CZ) gates at depth 235, vs 1658 CZ / depth 1188 for the dense/unrestricted
ansatz on the same backend. At ibm_kingston's own reported median CZ error
(0.20%), expected circuit fidelity is ~50%, a physically meaningful regime
(vs ~3% for the dense ansatz) -- this is why the hardware-native ansatz was
built at all rather than submitting the original circuit.
"""

import json
import time

import numpy as np
import qiskit
import ffsim
import ffsim.qiskit as ffq
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
from ffsim.qiskit.lucj_pass_manager import generate_lucj_pass_manager

from src.restricted_hamiltonian import load_restricted_hamiltonian
from src.lucj_ansatz_hw import LUCJ_HW_PATH

BACKEND_NAME = "ibm_kingston"
SHOTS = 10_000

JOB_META_PATH = "data/processed/hardware_job.json"
HARDWARE_COUNTS_PATH = "data/processed/hardware_counts.npz"


def load_hw_ansatz(path=LUCJ_HW_PATH):
    d = np.load(path)
    x_opt, n_reps, ncas = d["x_opt"], int(d["n_reps"]), int(d["ncas"])
    pairs_aa = [tuple(int(x) for x in p) for p in d["pairs_aa"]]
    pairs_ab = [tuple(int(x) for x in p) for p in d["pairs_ab"]]
    interaction_pairs = (pairs_aa, pairs_ab)
    op = ffsim.UCJOpSpinBalanced.from_parameters(
        x_opt, norb=ncas, n_reps=n_reps, interaction_pairs=interaction_pairs)
    return op, interaction_pairs, dict(e_hf=float(d["e_hf"]), e_opt=float(d["e_opt"]))


def build_logical_circuit(op, ncas, nelecas):
    qubits = qiskit.QuantumRegister(2 * ncas, name="q")
    circuit = qiskit.QuantumCircuit(qubits)
    circuit.append(ffq.PrepareHartreeFockJW(ncas, nelecas), qubits)
    circuit.append(ffq.UCJOpSpinBalancedJW(op), qubits)
    circuit.measure_all()
    return circuit


def build_isa_circuit(circuit, ncas, interaction_pairs, backend_name=BACKEND_NAME):
    service = QiskitRuntimeService()
    backend = service.backend(backend_name)
    pairs_aa, pairs_ab = interaction_pairs
    pm, accommodated_ab = generate_lucj_pass_manager(
        backend=backend, norb=ncas, connectivity="heavy-hex",
        interaction_pairs=(pairs_aa, pairs_ab), optimization_level=3,
    )
    isa_circuit = pm.run(circuit)
    return service, backend, isa_circuit, accommodated_ab


def submit_job(shots=SHOTS, backend_name=BACKEND_NAME):
    d_ham = load_restricted_hamiltonian()
    ncas, nelecas = d_ham["ncas"], d_ham["nelecas"]

    op, interaction_pairs, hw_info = load_hw_ansatz()
    circuit = build_logical_circuit(op, ncas, nelecas)
    service, backend, isa_circuit, accommodated_ab = build_isa_circuit(
        circuit, ncas, interaction_pairs, backend_name=backend_name)

    ops = isa_circuit.count_ops()
    n_2q = sum(v for k, v in ops.items() if k in ("cz", "ecr", "cx", "rzz"))
    print(f"backend={backend_name}, ISA circuit: depth={isa_circuit.depth()}, "
          f"size={isa_circuit.size()}, 2-qubit gates={n_2q}")
    print(f"pending jobs on {backend_name}: {backend.status().pending_jobs}")

    sampler = SamplerV2(mode=backend)
    job = sampler.run([isa_circuit], shots=shots)
    print(f"submitted job_id={job.job_id()} shots={shots}")

    meta = dict(job_id=job.job_id(), backend=backend_name, shots=shots,
                ncas=int(ncas), nelecas=list(int(x) for x in nelecas),
                n_2q_gates=int(n_2q), depth=int(isa_circuit.depth()),
                submitted_at=time.time(), e_hf=hw_info["e_hf"], e_opt=hw_info["e_opt"])
    with open(JOB_META_PATH, "w") as f:
        json.dump(meta, f, indent=2)
    print(f"saved job metadata -> {JOB_META_PATH}")
    return job, meta


def fetch_result(job_id=None, backend_name=BACKEND_NAME):
    if job_id is None:
        with open(JOB_META_PATH) as f:
            job_id = json.load(f)["job_id"]
    service = QiskitRuntimeService()
    job = service.job(job_id)
    print(f"job {job_id} status: {job.status()}")
    result = job.result()
    pub_result = result[0]
    bit_array = pub_result.data.meas
    return bit_array, job.status()


def load_hardware_bit_array(path=HARDWARE_COUNTS_PATH):
    """Reconstruct the exact BitArray qiskit-addon-sqd expects, from the raw
    packed-bit array saved by `fetch_result` -- NOT from the lossy
    bitstring/shots dict (`get_counts()`), so the SQD run downstream is fed
    the identical object type/packing as the simulator's `sample_lucj_circuit`
    output in src.sqd_benchmark.
    """
    from qiskit.primitives.containers.bit_array import BitArray
    d = np.load(path)
    return BitArray(d["packed_array"], int(d["num_bits"]))


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1 and sys.argv[1] == "fetch":
        bit_array, status = fetch_result()
        counts = bit_array.get_counts()
        print(f"status={status}, {len(counts)} unique bitstrings / "
              f"{bit_array.num_shots} shots")
        np.savez(HARDWARE_COUNTS_PATH,
                 packed_array=bit_array.array, num_bits=bit_array.num_bits,
                 bitstrings=np.array(list(counts.keys())),
                 shot_counts=np.array(list(counts.values())))
        print(f"saved -> {HARDWARE_COUNTS_PATH}")
    else:
        submit_job()

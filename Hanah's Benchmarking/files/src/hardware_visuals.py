"""One-time enrichment: real chip-topology and layout artifacts for the
`ibm_kingston` hardware run, generated while a live `QiskitRuntimeService()`
session is available so the output notebook itself never needs live IBM
credentials to render them.

Produces
--------
- `data/processed/hardware_error_map.png`: `qiskit.visualization.
  plot_error_map` snapshot of the full ibm_kingston chip (156 qubits) at
  the time of the run -- per-qubit readout error and per-edge two-qubit
  gate error, spatially laid out on the real heavy-hex topology. This is
  the actual noise landscape the job in src.hardware_run ran on, not a
  schematic.
- `data/processed/hardware_gate_map.png`: the same heavy-hex coupling
  graph with the physical qubits this specific job's ISA circuit was
  routed onto highlighted, so the two images together show both "what the
  whole chip looks like" and "which part of it we actually used".
- `data/processed/hardware_layout.npz`: the physical qubit indices used by
  the transpiled circuit (from its `layout.final_index_layout()`), plus
  the two-qubit-gate count/depth already saved by src.hardware_run, for
  the notebook to report without needing the live circuit object.
"""

import numpy as np
import matplotlib.pyplot as plt
from qiskit.visualization import plot_error_map, plot_gate_map

from src.restricted_hamiltonian import load_restricted_hamiltonian
from src.hardware_run import (
    load_hw_ansatz, build_logical_circuit, build_isa_circuit, BACKEND_NAME,
)

ERROR_MAP_PATH = "data/processed/hardware_error_map.png"
GATE_MAP_PATH = "data/processed/hardware_gate_map.png"
LAYOUT_PATH = "data/processed/hardware_layout.npz"


if __name__ == "__main__":
    d_ham = load_restricted_hamiltonian()
    ncas, nelecas = d_ham["ncas"], d_ham["nelecas"]

    op, interaction_pairs, _ = load_hw_ansatz()
    circuit = build_logical_circuit(op, ncas, nelecas)
    service, backend, isa_circuit, accommodated_ab = build_isa_circuit(
        circuit, ncas, interaction_pairs, backend_name=BACKEND_NAME)

    used_qubits = sorted(isa_circuit.layout.final_index_layout())
    print(f"physical qubits used by this job's ISA circuit "
          f"({len(used_qubits)} of {backend.num_qubits}): {used_qubits}")

    fig1 = plot_error_map(backend, figsize=(11, 9))
    fig1.savefig(ERROR_MAP_PATH, dpi=130, bbox_inches="tight")
    plt.close(fig1)
    print(f"saved -> {ERROR_MAP_PATH}")

    qubit_colors = ["#d62728" if q in used_qubits else "#cccccc"
                    for q in range(backend.num_qubits)]
    fig2 = plot_gate_map(backend, figsize=(11, 9), qubit_color=qubit_colors,
                          font_color="white")
    fig2.savefig(GATE_MAP_PATH, dpi=130, bbox_inches="tight")
    plt.close(fig2)
    print(f"saved -> {GATE_MAP_PATH}")

    props = backend.properties()
    cz_errors = np.array([
        param.value for gate in props.gates if gate.gate == "cz"
        for param in gate.parameters if param.name == "gate_error"
    ])
    median_cz_error = float(np.median(cz_errors))
    n2q = 338  # from src.hardware_transpile_comparison, hw-native ansatz on this backend
    expected_fidelity = (1 - median_cz_error) ** n2q
    print(f"median CZ error on {BACKEND_NAME}: {median_cz_error:.4f}, "
          f"expected fidelity at {n2q} 2q gates: {expected_fidelity:.3f}")

    np.savez(LAYOUT_PATH, used_qubits=np.array(used_qubits),
             num_qubits_backend=backend.num_qubits,
             median_cz_error=median_cz_error, expected_fidelity=expected_fidelity)
    print(f"saved -> {LAYOUT_PATH}")

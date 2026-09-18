"""Persist the logical-circuit diagram used in the notebook's real-hardware
section (data/processed/hardware_logical_circuit.png), which was previously
drawn ad hoc in an interactive shell. This script makes that figure
reproducible from a file under src/, consistent with every other figure in
this project.

Draws the LOGICAL circuit (Hartree-Fock state preparation + the
hardware-native LUCJ ansatz as single high-level ffsim gates + measurement),
not the ISA-transpiled circuit -- the transpiled depth/2-qubit-gate-count
numbers are reported separately by src.hardware_transpile_comparison.
"""

import matplotlib.pyplot as plt

from src.restricted_hamiltonian import load_restricted_hamiltonian
from src.hardware_run import load_hw_ansatz, build_logical_circuit

OUT_PATH = "data/processed/hardware_logical_circuit.png"


if __name__ == "__main__":
    d_ham = load_restricted_hamiltonian()
    ncas, nelecas = d_ham["ncas"], d_ham["nelecas"]

    op, interaction_pairs, _ = load_hw_ansatz()
    circuit = build_logical_circuit(op, ncas, nelecas)

    fig = circuit.draw("mpl", style="iqp", fold=-1)
    fig.savefig(OUT_PATH, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"saved -> {OUT_PATH}")

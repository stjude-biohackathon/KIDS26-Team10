"""Persist the side-by-side transpilation comparison that motivated
building a second, hardware-native LUCJ ansatz (src.lucj_ansatz_hw):
the ORIGINAL dense/unrestricted ansatz (src.lucj_ansatz) transpiled for
`ibm_kingston` with a generic preset pass manager, vs. the hardware-native
ansatz transpiled with ffsim's LUCJ-aware pass manager
(`generate_lucj_pass_manager`). Both numbers were first found by ad hoc
interactive exploration; this script makes them reproducible and gives
the notebook a file to load instead of a number asserted from memory.
"""

import numpy as np
import qiskit
import ffsim
import ffsim.qiskit as ffq
from qiskit_ibm_runtime import QiskitRuntimeService
from qiskit.transpiler import generate_preset_pass_manager
from ffsim.qiskit.lucj_pass_manager import generate_lucj_pass_manager

from src.restricted_hamiltonian import load_restricted_hamiltonian
from src.lucj_ansatz import LUCJ_PATH
from src.hardware_run import load_hw_ansatz, build_logical_circuit, BACKEND_NAME

OUT_PATH = "data/processed/hardware_transpile_comparison.npz"


def count_2q(ops):
    return sum(v for k, v in ops.items() if k in ("cz", "ecr", "cx", "rzz"))


if __name__ == "__main__":
    d_ham = load_restricted_hamiltonian()
    ncas, nelecas = d_ham["ncas"], d_ham["nelecas"]

    service = QiskitRuntimeService()
    backend = service.backend(BACKEND_NAME)

    # dense/unrestricted ansatz, generic transpilation
    d_dense = np.load(LUCJ_PATH)
    op_dense = ffsim.UCJOpSpinBalanced.from_parameters(
        d_dense["x_opt"], norb=ncas, n_reps=int(d_dense["n_reps"]))
    circuit_dense = build_logical_circuit(op_dense, ncas, nelecas)
    pm_generic = generate_preset_pass_manager(backend=backend, optimization_level=3)
    isa_dense = pm_generic.run(circuit_dense)
    ops_dense = isa_dense.count_ops()

    # hardware-native ansatz, LUCJ-aware transpilation
    op_hw, interaction_pairs, _ = load_hw_ansatz()
    circuit_hw = build_logical_circuit(op_hw, ncas, nelecas)
    pairs_aa, pairs_ab = interaction_pairs
    pm_lucj, _ = generate_lucj_pass_manager(
        backend=backend, norb=ncas, connectivity="heavy-hex",
        interaction_pairs=(pairs_aa, pairs_ab), optimization_level=3)
    isa_hw = pm_lucj.run(circuit_hw)
    ops_hw = isa_hw.count_ops()

    print(f"dense ansatz, generic PM on {BACKEND_NAME}: "
          f"depth={isa_dense.depth()}, 2q gates={count_2q(ops_dense)}")
    print(f"hw-native ansatz, LUCJ-aware PM on {BACKEND_NAME}: "
          f"depth={isa_hw.depth()}, 2q gates={count_2q(ops_hw)}")

    np.savez(
        OUT_PATH,
        depth_dense=isa_dense.depth(), n2q_dense=count_2q(ops_dense),
        depth_hw=isa_hw.depth(), n2q_hw=count_2q(ops_hw),
        backend=BACKEND_NAME,
    )
    print(f"saved -> {OUT_PATH}")

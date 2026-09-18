"""Hardware-native ("Local" in LUCJ) ansatz for the frozen, restricted
(10e,10o) active-space Hamiltonian: a second LUCJ variant, alongside the
dense one in src.lucj_ansatz, built to actually run on real superconducting
hardware rather than only in simulation.

Why a second ansatz is needed
------------------------------
src.lucj_ansatz's operator uses the DEFAULT (dense, all-to-all) diagonal
Coulomb connectivity: every pair of the 10 active orbitals interacts, in
both the same-spin and opposite-spin channels. Transpiling that circuit for
a real heavy-hex backend (ibm_kingston, verified directly: see the
docstring numbers below) requires so much SWAP-routing to realize the
implied all-to-all qubit connectivity that it explodes to depth 1188 and
1658 two-qubit (CZ) gates -- at ibm_kingston's own reported median CZ error
(0.20%, verified via backend.properties()), that is an expected circuit
fidelity of only ~3% at best (worse in practice, since routing will also
touch some of the backend's higher-error/dead qubits, max reported CZ
error = 100%). A bitstring distribution sampled from that circuit would be
dominated by noise, not by the ansatz.

This module builds the ansatz the way the LUCJ paper (Motta, Sung, Whaley,
Head-Gordon & Shee 2023, https://pubs.rsc.org/en/content/articlehtml/2023/
sc/d3sc02516k) actually intends for hardware: restrict the diagonal
Coulomb layer's `interaction_pairs` to a LINEAR CHAIN over the 10 orbitals
for the same-spin (alpha-alpha / beta-beta) sector, plus a small number of
opposite-spin (alpha-beta) connections matching ffsim's own documented
heavy-hex default (`[(p, p) for p in range(norb) if p % 4 == 0]`). ffsim's
gate-building code (`ffsim.qiskit.gates.diag_coulomb`) skips every zeroed
matrix entry when emitting gates (`if mat[i, j]: yield CircuitInstruction
(...)`, verified by reading its source directly), so restricting
interaction_pairs at the OPERATOR level directly removes those gates from
the circuit -- it is not just a hint to the transpiler. The orbital
rotation layer needs no such restriction: ffsim already decomposes any
orbital rotation into a linear-depth, nearest-neighbor-only Givens network
(Kivlichan et al. 2018), so it was never the source of the blowup.

ffsim additionally ships a purpose-built pass manager for exactly this
ansatz/topology pairing, `ffsim.qiskit.generate_lucj_pass_manager`, which
lays the two linear (alpha, beta) qubit chains onto a real heavy-hex
coupling graph via subgraph isomorphism instead of generic SWAP insertion.
Both pieces (restricted interaction_pairs + the LUCJ-aware pass manager)
are used together here and in the hardware submission script
(src.hardware_run) -- using only one of the two would still leave the
other bottleneck in place.

Same classical-optimization caveat as src.lucj_ansatz applies: no
converged single-reference (RHF/CCSD) amplitude source exists for this
active space (genuine open-shell singlet / diradical, see
src.spin_ladder), so this ansatz is also randomly initialized and
classically VQE-optimized from Hartree-Fock via ffsim's linear-method
optimizer -- still pure classical simulation, no quantum hardware
involved at this stage.
"""

import time

import numpy as np
import ffsim

from src.restricted_hamiltonian import load_restricted_hamiltonian
from src.lucj_ansatz import build_hamiltonian_operator

N_REPS = 1
MAX_ITER = 10
SEED = 0

LUCJ_HW_PATH = "data/processed/lucj_ansatz_hw.npz"


def linear_chain_interaction_pairs(norb):
    """Alpha-alpha/beta-beta: nearest-neighbor linear chain (norb-1 pairs).
    Alpha-beta: ffsim's own documented heavy-hex default, every 4th orbital
    (matches the ancilla-insertion pattern generate_lucj_pass_manager uses
    to route alpha-beta interactions on a heavy-hex coupling graph).
    """
    pairs_aa = [(p, p + 1) for p in range(norb - 1)]
    pairs_ab = [(p, p) for p in range(norb) if p % 4 == 0]
    return pairs_aa, pairs_ab


def optimize_lucj_hw(h1e, eri, ecore, ncas, nelecas, n_reps=N_REPS, max_iter=MAX_ITER,
                      seed=SEED):
    interaction_pairs = linear_chain_interaction_pairs(ncas)
    ham_op = build_hamiltonian_operator(h1e, eri, ecore, ncas, nelecas)
    vec_hf = ffsim.hartree_fock_state(ncas, nelecas)
    e_hf = float(np.vdot(vec_hf, ham_op @ vec_hf).real)

    ucj0 = ffsim.random.random_ucj_op_spin_balanced(
        ncas, n_reps=n_reps, interaction_pairs=interaction_pairs, seed=seed)
    x0 = ucj0.to_parameters(interaction_pairs=interaction_pairs)
    n_params = ffsim.UCJOpSpinBalanced.n_params(
        ncas, n_reps, interaction_pairs=interaction_pairs)
    assert x0.shape == (n_params,), (x0.shape, n_params)

    def params_to_vec(x):
        op = ffsim.UCJOpSpinBalanced.from_parameters(
            x, norb=ncas, n_reps=n_reps, interaction_pairs=interaction_pairs)
        return ffsim.apply_unitary(vec_hf, op, norb=ncas, nelec=nelecas)

    trace = []

    def cb(res):
        trace.append(float(res.fun))
        print(f"  iter {len(trace)}: E = {res.fun:.6f} Ha  ({time.time()-t0:.0f} s elapsed)",
              flush=True)

    t0 = time.time()
    e_random = float(np.vdot(params_to_vec(x0), ham_op @ params_to_vec(x0)).real)
    print(f"restricted interaction_pairs: aa/bb (linear chain) = {interaction_pairs[0]}")
    print(f"                              ab (heavy-hex default) = {interaction_pairs[1]}")
    print(f"n_params = {n_params} (vs {ffsim.UCJOpSpinBalanced.n_params(ncas, n_reps)} "
          f"for the unrestricted/dense ansatz)")
    print(f"E(HF determinant) = {e_hf:.6f} Ha")
    print(f"E(random hw-LUCJ, n_reps={n_reps}, seed={seed}) = {e_random:.6f} Ha")
    print(f"optimizing with linear-method VQE, max_iter={max_iter} ...")

    result = ffsim.optimize.minimize_linear_method(
        params_to_vec, ham_op, x0, maxiter=max_iter, callback=cb)
    dt = time.time() - t0

    x_opt = result.x
    op_opt = ffsim.UCJOpSpinBalanced.from_parameters(
        x_opt, norb=ncas, n_reps=n_reps, interaction_pairs=interaction_pairs)
    return dict(x_opt=x_opt, e_opt=float(result.fun), e_hf=e_hf, e_random=e_random,
                trace=np.array(trace), n_reps=n_reps, dt=dt, op=op_opt,
                interaction_pairs=interaction_pairs)


if __name__ == "__main__":
    d = load_restricted_hamiltonian()
    h1e, eri, ecore = d["h1e"], d["eri"], d["ecore"]
    ncas, nelecas = d["ncas"], d["nelecas"]
    e_casci = d["e_casci"]

    res = optimize_lucj_hw(h1e, eri, ecore, ncas, nelecas)

    print(f"\ndone in {res['dt']:.1f} s")
    print(f"E(HF)              = {res['e_hf']:.6f} Ha  "
          f"(gap to CASCI: {(res['e_hf']-e_casci)*1000:.2f} mHa)")
    print(f"E(random hw-LUCJ)  = {res['e_random']:.6f} Ha  "
          f"(gap to CASCI: {(res['e_random']-e_casci)*1000:.2f} mHa)")
    print(f"E(optimized hw-LUCJ) = {res['e_opt']:.6f} Ha  "
          f"(gap to CASCI: {(res['e_opt']-e_casci)*1000:.2f} mHa)")
    print(f"E_CASCI(S=0, exact)= {e_casci:.6f} Ha")
    print("optimization trace (Ha): " + ", ".join(f"{e:.4f}" for e in res["trace"]))

    pairs_aa, pairs_ab = res["interaction_pairs"]
    np.savez(
        LUCJ_HW_PATH,
        x_opt=res["x_opt"], n_reps=res["n_reps"],
        e_hf=res["e_hf"], e_random=res["e_random"], e_opt=res["e_opt"],
        trace=res["trace"], ncas=ncas, nelecas=np.array(nelecas),
        pairs_aa=np.array(pairs_aa), pairs_ab=np.array(pairs_ab),
    )
    print(f"\nsaved -> {LUCJ_HW_PATH}")

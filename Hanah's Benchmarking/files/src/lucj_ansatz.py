"""LUCJ (Local Unitary Cluster Jastrow) ansatz for the frozen, restricted
(10e,10o) active-space Hamiltonian, classically pre-optimized to give SQD a
non-trivial sampling distribution to work from.

Why not CCSD-seeded, as is standard practice
----------------------------------------------
The usual LUCJ recipe seeds `UCJOpSpinBalanced.from_t_amplitudes` with CCSD
t2 amplitudes from a converged RHF reference. Tried here first and it fails:
`pyscf.tools.fcidump.to_scf(...)` (via `ffsim.MolecularData.scf`) does not
converge in this active space -- checked explicitly, energy oscillates
wildly between SCF cycles (delta_E ~ 0.1-1 Ha at cycle 50, never settling).
This is the expected signature of a genuine open-shell singlet / diradical
active space (the same physics that gives the S=0..5 spin ladder in
src.spin_ladder its sub-mHa splittings): no single closed-shell determinant
is a good zeroth-order reference, so single-reference RHF/CCSD is not a
valid amplitude source here.

What is done instead
---------------------
The LUCJ ansatz is built with randomly initialized parameters (fixed seed,
reproducible) and then classically variationally optimized against the
*exact* frozen active-space Hamiltonian using ffsim's linear-method VQE
optimizer (`ffsim.optimize.minimize_linear_method`, Motta et al. 2014),
starting from the Hartree-Fock (aufbau) reference determinant. This is
still classical simulation throughout (statevector, via ffsim's fast
Jordan-Wigner-free simulator) -- there is no quantum hardware or quantum
resource involved at this stage, exactly as in a normal VQE pre-training
step before using an ansatz to generate hardware/sampler counts.

Compute budget: each linear-method iteration costs ~30-40 s on this
2-core sandbox for n_reps=1 (210 parameters), so MAX_ITER is deliberately
kept small; the resulting ansatz is under-converged relative to a full VQE
optimum (see the printed energy trace) but is NOT required to be
variationally exact for the SQD benchmark: SQD's self-consistent
configuration recovery (src.sqd_benchmark) is specifically designed to
extract a much better energy than the ansatz's own energy, as long as the
sampled bitstring distribution has non-negligible overlap with the true
ground-state determinants. That is the property being tested here, not
ansatz optimality.
"""

import time

import numpy as np
import ffsim

from src.restricted_hamiltonian import load_restricted_hamiltonian

N_REPS = 1
MAX_ITER = 10
SEED = 0

LUCJ_PATH = "data/processed/lucj_ansatz.npz"


def build_hamiltonian_operator(h1e, eri, ecore, ncas, nelecas):
    mol_ham = ffsim.MolecularHamiltonian(h1e, eri, constant=ecore)
    return ffsim.linear_operator(mol_ham, norb=ncas, nelec=nelecas)


def optimize_lucj(h1e, eri, ecore, ncas, nelecas, n_reps=N_REPS, max_iter=MAX_ITER, seed=SEED):
    ham_op = build_hamiltonian_operator(h1e, eri, ecore, ncas, nelecas)
    vec_hf = ffsim.hartree_fock_state(ncas, nelecas)
    e_hf = float(np.vdot(vec_hf, ham_op @ vec_hf).real)

    ucj0 = ffsim.random.random_ucj_op_spin_balanced(ncas, n_reps=n_reps, seed=seed)
    x0 = ucj0.to_parameters()

    def params_to_vec(x):
        op = ffsim.UCJOpSpinBalanced.from_parameters(x, norb=ncas, n_reps=n_reps)
        return ffsim.apply_unitary(vec_hf, op, norb=ncas, nelec=nelecas)

    trace = []

    def cb(res):
        trace.append(float(res.fun))
        print(f"  iter {len(trace)}: E = {res.fun:.6f} Ha  ({time.time()-t0:.0f} s elapsed)",
              flush=True)

    t0 = time.time()
    e_random = float(np.vdot(params_to_vec(x0), ham_op @ params_to_vec(x0)).real)
    print(f"E(HF determinant) = {e_hf:.6f} Ha")
    print(f"E(random LUCJ, n_reps={n_reps}, seed={seed}) = {e_random:.6f} Ha")
    print(f"optimizing with linear-method VQE, max_iter={max_iter} ...")

    result = ffsim.optimize.minimize_linear_method(
        params_to_vec, ham_op, x0, maxiter=max_iter, callback=cb)
    dt = time.time() - t0

    x_opt = result.x
    op_opt = ffsim.UCJOpSpinBalanced.from_parameters(x_opt, norb=ncas, n_reps=n_reps)
    return dict(x_opt=x_opt, e_opt=float(result.fun), e_hf=e_hf, e_random=e_random,
                trace=np.array(trace), n_reps=n_reps, dt=dt, op=op_opt)


if __name__ == "__main__":
    d = load_restricted_hamiltonian()
    h1e, eri, ecore = d["h1e"], d["eri"], d["ecore"]
    ncas, nelecas = d["ncas"], d["nelecas"]
    e_casci = d["e_casci"]

    res = optimize_lucj(h1e, eri, ecore, ncas, nelecas)

    print(f"\ndone in {res['dt']:.1f} s")
    print(f"E(HF)              = {res['e_hf']:.6f} Ha  "
          f"(gap to CASCI: {(res['e_hf']-e_casci)*1000:.2f} mHa)")
    print(f"E(random LUCJ)     = {res['e_random']:.6f} Ha  "
          f"(gap to CASCI: {(res['e_random']-e_casci)*1000:.2f} mHa)")
    print(f"E(optimized LUCJ)  = {res['e_opt']:.6f} Ha  "
          f"(gap to CASCI: {(res['e_opt']-e_casci)*1000:.2f} mHa)")
    print(f"E_CASCI(S=0, exact)= {e_casci:.6f} Ha")
    print("optimization trace (Ha): " + ", ".join(f"{e:.4f}" for e in res["trace"]))

    np.savez(
        LUCJ_PATH,
        x_opt=res["x_opt"], n_reps=res["n_reps"],
        e_hf=res["e_hf"], e_random=res["e_random"], e_opt=res["e_opt"],
        trace=res["trace"], ncas=ncas, nelecas=np.array(nelecas),
    )
    print(f"\nsaved -> {LUCJ_PATH}")

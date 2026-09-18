"""Frozen, SPIN-RESTRICTED (10e,10o) active-space Hamiltonian for the
SQD/LUCJ benchmark and the matched-dimension sampler baseline.

Why this is a SEPARATE Hamiltonian from src.hamiltonian
--------------------------------------------------------
src.hamiltonian freezes the UNRESTRICTED (BS-localized) active-space
Hamiltonian needed for the broken-symmetry picture (src.spin_observables).
qiskit-addon-sqd and ffsim's LUCJ ansatz (UCJOpSpinBalanced,
ffsim.MolecularData) are built for a single common one-/two-body integral
set shared by both spins (`one_body_integrals`, `two_body_integrals`: see
ffsim.MolecularData -- no separate alpha/beta arrays), i.e. they assume a
spin-restricted active space. Feeding them the BS Hamiltonian's h1e_a != h1e_b
integrals is not an option the API supports, and conceptually SQD papers
benchmark against ordinary (spin-restricted) CASCI, not broken-symmetry
CASCI.

We reuse the SAME common (spin-symmetric) active orbitals already validated
in src.spin_ladder (from AVAS on the HIGH-SPIN UHF reference, which preserves
Fe1<->Fe2 symmetry -- see that module's docstring for why). The restricted
CASCI(10,(5,5)) ground state in this active space is the S=0 state of the
spin ladder already computed there (E=-5013.45063122 Ha, <S^2>=0 exactly):
this module just re-derives and freezes the same physics as one- and
two-body integrals for ffsim/SQD to consume, and re-verifies the energy
against that already-computed number.
"""

import time

import numpy as np
from pyscf import ao2mo, mcscf

from src.meanfield import build_mol
from src.spin_ladder import build_common_active_orbitals, compute_spin_ladder, NCAS, NELECAS

FROZEN_RESTRICTED_HAM_PATH = "data/processed/restricted_active_space_hamiltonian.npz"


def build_and_freeze_restricted_hamiltonian(mol, mo_common, ncas=NCAS, nelecas=NELECAS,
                                             out_path=FROZEN_RESTRICTED_HAM_PATH):
    """Freeze h1e/eri/ecore for the restricted (10,10) active space, and
    separately pin down the true S=0 ground-state energy via
    src.spin_ladder.compute_spin_ladder's tight-convergence diagonalization.

    NOTE: a plain `mcscf.CASCI(...).kernel()` with default (loose) FCI
    convergence on this system lands on a spin-CONTAMINATED state (checked:
    S^2=3.615, not any clean S(S+1)) 0.22 mHa above the true S=0 ground
    state, for the same near-degenerate-spin-ladder reason documented in
    src.spin_ladder -- so it must not be used as the SQD benchmark target.
    get_h1eff/get_h2eff only depend on the *orbitals* (mo_common), not on
    which CI root was found, so they are unaffected by this and are computed
    directly with a throwaway CASCI object below.
    """
    mc = mcscf.CASCI(mol, ncas, nelecas)
    ncore = mc.ncore  # inferred from (mol.nelectron, mol.spin, nelecas); no kernel() needed
    mo_cas = mo_common[:, ncore:ncore + ncas]

    h1e, ecore = mc.get_h1eff(mo_common)
    eri = ao2mo.restore(1, mc.get_h2eff(mo_common), ncas)

    ladder, mc_ladder = compute_spin_ladder(mol, mo_common, ncas=ncas, nelecas=nelecas)
    e_s0, ss_s0 = ladder[0]
    assert abs(ss_s0) < 1e-6, f"S=0 state is not spin-pure: <S^2>={ss_s0}"

    np.savez(
        out_path,
        h1e=h1e, eri=eri, ecore=ecore,
        mo_coeff=mo_common, mo_cas=mo_cas, ncore=ncore,
        ncas=ncas, nelecas=np.array(nelecas),
        e_casci=e_s0,
    )
    return mc_ladder, h1e, eri, ecore, e_s0


def load_restricted_hamiltonian(path=FROZEN_RESTRICTED_HAM_PATH):
    d = np.load(path)
    return dict(
        h1e=d["h1e"], eri=d["eri"], ecore=float(d["ecore"]),
        mo_coeff=d["mo_coeff"], mo_cas=d["mo_cas"], ncore=int(d["ncore"]),
        ncas=int(d["ncas"]), nelecas=tuple(int(x) for x in d["nelecas"]),
        e_casci=float(d["e_casci"]),
    )


if __name__ == "__main__":
    mol, _ = build_mol()

    t0 = time.time()
    mo_common, ncas, mf_hs = build_common_active_orbitals()
    print(f"common active orbitals: ncas={ncas}  ({time.time()-t0:.1f} s)")

    t0 = time.time()
    mc, h1e, eri, ecore, e_s0 = build_and_freeze_restricted_hamiltonian(mol, mo_common)
    dt = time.time() - t0
    print(f"restricted CASCI({NCAS}e,{NCAS}o) built and frozen ({dt:.1f} s)")
    print(f"  E_CASCI(S=0) = {e_s0:.8f} Ha")
    print("  verification: matches src.spin_ladder S=0 rung E=-5013.45063122 Ha "
          f"to {abs(e_s0 - (-5013.45063122)):.2e} Ha")
    print(f"  h1e shape {h1e.shape}, eri shape {eri.shape}, ecore={ecore:.8f} Ha")
    print(f"saved -> {FROZEN_RESTRICTED_HAM_PATH}")

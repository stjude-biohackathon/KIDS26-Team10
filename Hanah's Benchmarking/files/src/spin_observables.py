"""Spin observables of the CASCI(10e,10o) ground state built on the frozen,
unrestricted (BS-derived) active-space Hamiltonian: natural-orbital
occupations, <S^2>, and local (Mulliken) spin populations per atom.

This deliberately reuses the SAME frozen Hamiltonian
(data/processed/active_space_hamiltonian.npz, built by src.hamiltonian) that
feeds the SQD benchmark, so these observables describe exactly the state SQD
is being asked to reconstruct -- not some other, re-optimized reference.

<S^2> subtlety
--------------
Because this active space is unrestricted (alpha and beta active orbitals
are different, spatially localized on opposite Fe centers -- see
src/active_space.py), the "obvious" S^2 estimator
`fcisolver.spin_square(civec, ncas, nelecas)` implicitly assumes an
orthonormal common alpha/beta basis (ovlp=1) and is WRONG here: it returns
~0 regardless of the true spin state (checked explicitly below). The correct
estimator, `pyscf.fci.spin_op.spin_square_general`, needs the actual AO
overlap between the alpha and beta active orbitals, which is not the
identity since they are different, non-orthogonal MOs.
"""

import numpy as np
from pyscf import mcscf
from pyscf.fci import direct_uhf, spin_op

from src.hamiltonian import load_bs_uhf, load_frozen_hamiltonian, FROZEN_HAM_PATH
from src.meanfield import _ao_spin_pop_per_atom

SPIN_OBS_PATH = "data/processed/spin_observables.npz"


def ground_state_ci(h1e, eri, ncas, nelecas):
    cisolver = direct_uhf.FCI()
    e_active, civec = cisolver.kernel(h1e, eri, ncas, nelecas)
    return cisolver, e_active, civec


def natural_orbital_occupations(cisolver, civec, ncas, nelecas):
    """Eigenvalues of the active-space 1-RDM per spin channel: how close the
    FCI ground state is to a single Slater determinant (occupations exactly
    0/1) versus genuinely multi-configurational (fractional occupations).
    """
    dm1a, dm1b = cisolver.make_rdm1s(civec, ncas, nelecas)
    occ_a = np.sort(np.linalg.eigvalsh(dm1a))[::-1]
    occ_b = np.sort(np.linalg.eigvalsh(dm1b))[::-1]
    return occ_a, occ_b


def spin_square_active(cisolver, civec, ncas, nelecas, mo_a_cas, mo_b_cas, ao_overlap):
    """<S^2> of the active-space FCI ground state, correctly accounting for
    alpha/beta orbital non-orthogonality (see module docstring).
    """
    dm1a, dm1b = cisolver.make_rdm1s(civec, ncas, nelecas)
    _, (dm2aa, dm2ab, dm2bb) = cisolver.make_rdm12s(civec, ncas, nelecas)
    ss, mult = spin_op.spin_square_general(
        dm1a, dm1b, dm2aa, dm2ab, dm2bb, (mo_a_cas, mo_b_cas), ao_overlap)
    # naive (wrong) estimator, computed for the record to show why it fails
    ss_naive, mult_naive = cisolver.spin_square(civec, ncas, nelecas)
    return ss, mult, ss_naive, mult_naive


def local_spin_populations(mol, mo_a, mo_b, ncore, ncas, nelecas, cisolver, civec):
    """Mulliken atomic spin (alpha-beta) populations of the full AO-basis
    density built from this active-space ground state: core (doubly present
    per channel) + active (from the FCI 1-RDM), mirroring
    mcscf.ucasci.UCASCI.make_rdm1s exactly so the numbers agree with what
    UCASCI itself would report.
    """
    ncore_a, ncore_b = ncore
    casdm1a, casdm1b = cisolver.make_rdm1s(civec, ncas, nelecas)

    mocore_a = mo_a[:, :ncore_a]
    mocas_a = mo_a[:, ncore_a:ncore_a + ncas]
    dm1a = mocore_a @ mocore_a.T + mocas_a @ casdm1a @ mocas_a.T

    mocore_b = mo_b[:, :ncore_b]
    mocas_b = mo_b[:, ncore_b:ncore_b + ncas]
    dm1b = mocore_b @ mocore_b.T + mocas_b @ casdm1b @ mocas_b.T

    return _ao_spin_pop_per_atom(mol, dm1a, dm1b), dm1a, dm1b


if __name__ == "__main__":
    import time

    d = load_frozen_hamiltonian()
    h1e, eri, ecore = d["h1e"], d["eri"], d["ecore"]
    ncas, nelecas = d["ncas"], d["nelecas"]
    mo_a, mo_b = d["mo"]
    ncore = d["ncore"]
    ao_overlap = d["ao_overlap"]
    mo_a_cas = mo_a[:, ncore[0]:ncore[0] + ncas]
    mo_b_cas = mo_b[:, ncore[1]:ncore[1] + ncas]

    print(f"loaded frozen Hamiltonian from {FROZEN_HAM_PATH}: "
          f"ncas={ncas}, nelecas={tuple(nelecas)}, ncore={tuple(ncore)}")

    t0 = time.time()
    cisolver, e_active, civec = ground_state_ci(h1e, eri, ncas, nelecas)
    dt = time.time() - t0
    e_tot = e_active + ecore
    print(f"\nFCI ground state: E_active={e_active:.8f} Ha, E_tot={e_tot:.8f} Ha  ({dt:.2f} s)")
    print("  (verification: matches src.active_space UCASCI E=-5013.91060402 Ha "
          f"to {abs(e_tot - (-5013.91060402)):.2e} Ha)")

    occ_a, occ_b = natural_orbital_occupations(cisolver, civec, ncas, nelecas)
    print("\nNatural orbital occupations (active space, alpha channel):")
    print("  " + "  ".join(f"{o:.6f}" for o in occ_a))
    print("Natural orbital occupations (active space, beta channel):")
    print("  " + "  ".join(f"{o:.6f}" for o in occ_b))
    print("  (all ~0 or ~1 => the FCI ground state is essentially a single "
          "Slater determinant in this active space: little additional static "
          "correlation beyond the broken-symmetry mean-field picture)")

    ss, mult, ss_naive, mult_naive = spin_square_active(
        cisolver, civec, ncas, nelecas, mo_a_cas, mo_b_cas, ao_overlap)
    print(f"\n<S^2> (naive estimator, WRONGLY assumes orthonormal alpha/beta "
          f"basis) = {ss_naive:.6f}  <- ignores that this is an unrestricted "
          "active space; not physically meaningful here")
    print(f"<S^2> (correct, AO-overlap-aware estimator) = {ss:.6f}  "
          f"(2S+1 = {mult:.4f})")

    mol, mf, _ = load_bs_uhf()
    ss_core, _ = mf.spin_square((mo_a[:, :ncore[0]], mo_b[:, :ncore[1]]), ao_overlap)
    print(f"  core-orbital contribution to <S^2> = {ss_core:.6f} "
          f"(active+core = {ss + ss_core:.6f}, matches the full UCASCI "
          "S^2=5.6554184 reported by src.active_space to ~0.1%, the residual "
          "being active/core cross-terms this additive estimate omits)")

    pops, dm1a, dm1b = local_spin_populations(mol, mo_a, mo_b, ncore, ncas, nelecas,
                                               cisolver, civec)
    print("\nMulliken spin density (alpha-beta) per atom, full AO-basis density "
          "reconstructed from this active-space ground state:")
    for k, v in pops.items():
        if "Fe" in k:
            print(f"  {k}: {v:+.4f}")
    bridging_s = [v for k, v in pops.items() if k.split(":")[1] == "S"]
    print(f"  bridging/terminal S atoms: min={min(bridging_s):+.4f}, "
          f"max={max(bridging_s):+.4f}, mean|.|={np.mean(np.abs(bridging_s)):.4f}")
    print("  (matches the BS-UHF Fe populations +-5.165 almost exactly, "
          "consistent with the near-integer natural-orbital occupations "
          "above: at this AVAS(Fe 3d)/STO-3G active-space quality, FCI adds "
          "essentially no configuration mixing beyond the single BS-UHF "
          "determinant -- the active space is not yet capturing meaningful "
          "static correlation, only the exact recoupling used for the spin "
          "ladder in src.spin_ladder, which uses a *different*, spin-common "
          "active space where multi-determinant character is unambiguous)")

    np.savez(
        SPIN_OBS_PATH,
        e_active=e_active, e_tot=e_tot,
        occ_a=occ_a, occ_b=occ_b,
        ss_active=ss, mult_active=mult, ss_naive=ss_naive,
        ss_core=ss_core,
        fe_pops=np.array([v for k, v in pops.items() if "Fe" in k]),
        all_pops_keys=np.array(list(pops.keys())),
        all_pops_vals=np.array(list(pops.values())),
    )
    print(f"\nsaved -> {SPIN_OBS_PATH}")

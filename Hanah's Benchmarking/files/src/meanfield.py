"""Broken-symmetry UHF/def2-SVP mean field for [Fe2S2(SMe)4]2-.

Electron-counting rationale (documented so the charge/spin choice is
checkable, not asserted):
    4 x SMe-  ligand charge  = -4
    2 x S(2-) bridging sulfido = -4
    total ligand charge = -8
    overall complex charge = -2  =>  Fe2 charge = -2 - (-8) = +6
    => both Fe formally Fe(III), d5 high-spin (S=5/2 each) in the localized
       picture, antiferromagnetically coupled to a total S=0 ground state.
This is exactly the oxidized [2Fe-2S]2+ core studied in Noodleman's original
broken-symmetry (BS) DFT papers (Noodleman, J. Chem. Phys. 1981, 74, 5737),
the standard benchmark system for BS mean-field theory of Fe-S clusters.

Mean-field recipe: UHF with charge=-2, spin=0 (Ms=0), started from a
symmetry-broken initial density (Fe1 excess alpha / Fe2 excess beta) so SCF
converges to the antiferromagnetic broken-symmetry solution rather than the
(unphysical, higher-energy) symmetric one.
"""

import numpy as np
from pyscf import gto, scf

from src.geometry import build_fe2s2_sme4, atoms_to_pyscf_atom

BASIS = "sto-3g"
# NOTE: STO-3G (122 basis functions here) is chosen over def2-SVP (286
# functions) purely for tractability: this sandbox has 2 CPU cores, and
# DF-UHF/def2-SVP on this 186-electron molecule did not converge in 25+
# minutes even for the easier high-spin reference. STO-3G is a genuine
# accuracy tradeoff (no polarization functions, minimal valence), so
# absolute energies here are a methodology demonstration of the
# AVAS -> CASCI -> LUCJ/SQD workflow, not publication-accuracy numbers.
# The active-space Hamiltonian this produces is frozen and reused
# identically by both the CASCI spin-ladder validation and the SQD
# benchmark, so any basis-set error affects both sides equally.
CHARGE = -2
SPIN = 0        # Ms = 0 target (Nalpha = Nbeta), the broken-symmetry singlet
SPIN_HS = 10    # Ms = 5 high-spin reference (2 x Fe(III) d5, S=5/2, aligned)


def build_mol(basis=BASIS, charge=CHARGE, spin=SPIN, verbose=3):
    atoms, geom_checks = build_fe2s2_sme4(charge=charge)
    mol = gto.M(
        atom=atoms_to_pyscf_atom(atoms),
        basis=basis,
        charge=charge,
        spin=spin,
        verbose=verbose,
    )
    return mol, geom_checks


def _fe_indices(mol):
    fe_idx = [i for i, (sym, _) in enumerate(mol._atom) if sym == "Fe"]
    assert len(fe_idx) == 2, f"expected 2 Fe atoms, found {len(fe_idx)}"
    return fe_idx


def _ao_spin_pop_per_atom(mol, dma, dmb):
    """Mulliken atomic spin populations (alpha - beta) from AO-basis density
    matrices: pop_mu = sum_nu (dm_spin @ S)_{mu mu}, summed over each atom's
    AOs.
    """
    dm_spin = dma - dmb
    S = mol.intor_symmetric("int1e_ovlp")
    ao_pop = np.einsum("ij,ji->i", dm_spin, S)
    per_atom = {}
    for i, (sym, _) in enumerate(mol._atom):
        _, _, ao_start, ao_end = mol.aoslice_by_atom()[i]
        per_atom[f"{i}:{sym}"] = float(ao_pop[ao_start:ao_end].sum())
    return per_atom


def _flip_fe2_density(mol, dma, dmb, fe2_ao):
    """Exchange the alpha/beta AO density on Fe2's diagonal (Fe2-Fe2) block
    only, the standard 'spin-flip' recipe for building a broken-symmetry
    initial guess from a stable reference density. Off-diagonal (Fe2 x rest)
    blocks are left untouched: swapping them too would require also
    reassigning the coupling terms, and the diagonal-block swap alone is
    enough to seed the opposite local moment on Fe2 for SCF to relax from.
    """
    idx = np.ix_(fe2_ao, fe2_ao)
    dma_new, dmb_new = dma.copy(), dmb.copy()
    dma_new[idx] = dmb[idx]
    dmb_new[idx] = dma[idx]
    return dma_new, dmb_new


def run_high_spin_uhf(mol_hs, max_cycle=100, conv_tol=1e-8, density_fit=True):
    """Converge the Ms=5 high-spin UHF reference (both Fe(III) d5 aligned).
    A single, well-separated Slater determinant: no near-degeneracy, so
    plain DIIS-SCF converges reliably. Used only to build a chemically
    sound initial density for the broken-symmetry Ms=0 calculation.
    """
    mf = scf.UHF(mol_hs)
    if density_fit:
        mf = mf.density_fit()
    mf.max_cycle = max_cycle
    mf.conv_tol = conv_tol
    mf.kernel()
    return mf


def run_bs_uhf(mol, max_cycle=100, conv_tol=1e-8, density_fit=True, use_newton=True):
    """Run Ms=0 UHF seeded from the Fe2-spin-flipped high-spin density, to
    land on the antiferromagnetic broken-symmetry (BS) solution rather than
    the symmetric/delocalized one.

    density_fit=True uses the resolution-of-identity (RI-JK) approximation
    for the two-electron integrals: on this machine's 2 CPU cores, exact
    4-index UHF/def2-SVP for this 286-basis-function, 186-electron molecule
    is impractically slow; RI-JK cuts the integral cost from O(N^4) to
    O(N^2 Naux) with a standard, well-tested auxiliary basis, at the cost of
    a controlled sub-mHa RI error. This only affects which orbitals
    AVAS/CASCI start from, not the final correlated energies.

    use_newton=True switches to PySCF's second-order (Newton-Raphson) SCF
    solver for the Ms=0 step: plain first-order DIIS-SCF tends to slide back
    to the higher-symmetry, lower-<S^2> stationary point even from a
    symmetry-broken guess; Newton-SCF converges to the stationary point
    nearest the starting density, which is what an initial guess is for.
    """
    fe_idx = _fe_indices(mol)

    mol_hs, _ = build_mol(basis=mol.basis, charge=mol.charge, spin=SPIN_HS, verbose=0)
    mf_hs = run_high_spin_uhf(mol_hs, density_fit=density_fit)

    fe2_ao = mol.search_ao_label(f"{fe_idx[1]} Fe")
    dma_hs, dmb_hs = mf_hs.make_rdm1()
    dma0, dmb0 = _flip_fe2_density(mol, dma_hs, dmb_hs, fe2_ao)

    # Verify the flip actually broke symmetry *before* spending an SCF on
    # it: Fe1 and Fe2 spin populations in the seed density must have
    # opposite sign, or this whole routine is pointless.
    seed_spins = _ao_spin_pop_per_atom(mol, dma0, dmb0)
    fe_spins = [v for k, v in seed_spins.items() if "Fe" in k]
    print(f"seed (pre-SCF) Fe spin populations: {[f'{v:+.3f}' for v in fe_spins]}")
    assert fe_spins[0] * fe_spins[1] < 0, (
        f"initial-guess density flip did not break symmetry: Fe spins {fe_spins} "
        "have the same sign; _flip_fe2_density is broken")

    mf = scf.UHF(mol)
    if density_fit:
        mf = mf.density_fit()
    mf.max_cycle = max_cycle
    mf.conv_tol = conv_tol

    if use_newton:
        mf = mf.newton()
        mf.max_cycle_inner = 100

    mf.chkfile = "data/processed/bs_uhf.chk"
    mf.kernel(dm0=(dma0, dmb0))
    return mf, fe_idx, mf_hs


def spin_density_per_atom(mf, mol):
    """Mulliken atomic spin populations (alpha - beta) of the converged
    state, for checking it is genuinely broken-symmetry (opposite sign on
    the two Fe centers).
    """
    dm_a, dm_b = mf.make_rdm1()
    return _ao_spin_pop_per_atom(mol, dm_a, dm_b)


if __name__ == "__main__":
    import time

    mol, geom_checks = build_mol()
    print(f"nao = {mol.nao}, nelectron = {mol.nelectron}, spin(2S) = {mol.spin}")

    t0 = time.time()
    mf, fe_idx, mf_hs = run_bs_uhf(mol)
    dt = time.time() - t0
    print(f"high-spin (Ms=5) UHF: E = {mf_hs.e_tot:.8f} Ha, converged={mf_hs.converged}")
    print(f"BS-UHF (Ms=0) converged: {mf.converged}   E(UHF) = {mf.e_tot:.8f} Ha   ({dt:.1f} s)")
    print(f"E(BS) - E(HS) = {(mf.e_tot - mf_hs.e_tot)*1000:.3f} mHa "
          "(expect < 0: antiferromagnetic coupling favors the low-spin state)")

    ss, s_mult = mf.spin_square()
    print(f"<S^2> = {ss:.4f}  (2S+1 = {s_mult:.4f}; ideal Ms=0 singlet has <S^2>=0, "
          "a genuine BS solution for 2x S=5/2 sites has <S^2> >> 0, roughly 5-15)")

    spins = spin_density_per_atom(mf, mol)
    print("Mulliken spin density (alpha-beta) on Fe atoms (expect opposite signs, |.|~3-4):")
    for k, v in spins.items():
        if "Fe" in k:
            print(f"  atom {k}: {v:+.3f}")

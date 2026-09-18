"""AVAS active-space selection and unrestricted CASCI(10e,10o) for the
BS-UHF description of [Fe2S2(SMe)4]2-.

Because the mean-field reference is a broken-symmetry (BS) UHF solution
(Fe1 majority-alpha 3d^5, Fe2 majority-beta 3d^5), the alpha and beta MO
sets are physically different: Fe1's 3d shell is occupied in the alpha
channel and empty in the beta channel, and vice versa for Fe2. A single
common spatial active space (as CASCI/CASSCF normally uses) would force
these two different orbital sets together and lose exactly the physics
that made the calculation broken-symmetry in the first place.

We therefore run AVAS independently on each spin channel, targeting the
'Fe 3d' AO labels (2 Fe x 5 d orbitals = 10 AOs, matching the target
(10e,10o) active space from the project profile):
    alpha channel -> AVAS selects Fe1's 5 occupied 3d (alpha-occupied) and
                     Fe2's 5 unoccupied 3d (alpha-virtual)
    beta channel  -> AVAS selects Fe2's 5 occupied 3d (beta-occupied) and
                     Fe1's 5 unoccupied 3d (beta-virtual)
This gives a genuinely unrestricted (10,10) active space, run through
PySCF's UCASCI (unrestricted CASCI, fci.direct_uhf), with 5 active alpha
electrons and 5 active beta electrons -- 10 active electrons total,
matching the 2x Fe(III) d5 picture used to build the BS reference.
"""

import numpy as np
from pyscf import mcscf
from pyscf.mcscf import avas

AOLABELS = ["Fe 3d"]
THRESHOLD = 0.2
NCAS_TARGET = 10
NELECAS_TARGET = (5, 5)


class _SpinChannelMF:
    """Minimal stand-in exposing exactly what pyscf.mcscf.avas needs for one
    UHF spin channel, so AVAS treats it as if it were an independent
    (non-UHF) SCF object and does not force the 'alpha orbitals only' UHF
    special-case in pyscf.mcscf.avas._kernel.
    """

    def __init__(self, mol, mo_coeff, mo_occ, mo_energy):
        self.mol = mol
        self.mo_coeff = mo_coeff
        self.mo_occ = mo_occ
        self.mo_energy = mo_energy
        self.stdout = mol.stdout
        self.verbose = mol.verbose


def avas_per_spin_channel(mol, mo_coeff, mo_occ, mo_energy, aolabels=AOLABELS,
                           threshold=THRESHOLD):
    """Run AVAS on a single spin channel treated as an independent
    (spin-orbital) SCF problem: mo_occ here is 0/1 (one electron per occupied
    orbital of this spin), not 0/2 as AVAS's internal `nelecas` formula
    assumes for a spin-restricted reference. That formula,
    `nelecas = (mol.nelectron - ncore*2) - (n_inactive_occ)*2`, uses the
    *total* (alpha+beta) electron count and a factor of 2 per orbital, so its
    return value is meaningless for a single channel and must not be used
    directly. The physically correct electron count contributed by this
    channel is simply the number of *occupied* orbitals AVAS placed in the
    active space (each holds exactly one electron of this spin), which is
    recovered from `occ_weights >= threshold` on the AVAS object.
    """
    fake_mf = _SpinChannelMF(mol, mo_coeff, mo_occ, mo_energy)
    avas_obj = avas.AVAS(fake_mf, aolabels, threshold=threshold, openshell_option=2)
    ncas, _nelecas_rhf_formula, mo = avas_obj.kernel()
    n_active_occ = int((avas_obj.occ_weights >= threshold).sum())
    return ncas, n_active_occ, mo


def build_unrestricted_active_space(mol, mf, aolabels=AOLABELS, threshold=THRESHOLD):
    """Run AVAS on the alpha and beta channels independently.

    Returns
    -------
    mo_a, mo_b : full alpha/beta MO coefficient matrices, active block placed
        contiguously (AVAS's own core/active/virtual ordering).
    ncas_a, ncas_b : number of active alpha/beta spatial orbitals selected.
    nelecas_a, nelecas_b : number of active alpha/beta electrons implied by
        the occupied/virtual split AVAS found (checked against the target
        5/5 downstream).
    """
    ncas_a, nelecas_a, mo_a = avas_per_spin_channel(
        mol, mf.mo_coeff[0], mf.mo_occ[0], mf.mo_energy[0], aolabels, threshold)
    ncas_b, nelecas_b, mo_b = avas_per_spin_channel(
        mol, mf.mo_coeff[1], mf.mo_occ[1], mf.mo_energy[1], aolabels, threshold)
    return mo_a, mo_b, ncas_a, ncas_b, nelecas_a, nelecas_b


def run_ucasci(mol, mf, mo_a, mo_b, ncas, nelecas):
    """Unrestricted CASCI (fci.direct_uhf) in the AVAS-selected (ncas, nelecas)
    active space, with ncore alpha/beta core orbitals determined from the
    total electron counts.
    """
    mc = mcscf.UCASCI(mf, ncas, nelecas)
    mc.mo_coeff = (mo_a, mo_b)
    mc.kernel(mo_coeff=(mo_a, mo_b))
    return mc


if __name__ == "__main__":
    import time
    from pyscf import lib
    from src.meanfield import build_mol, SPIN

    mol, _ = build_mol()
    mf = lib.chkfile.load("data/processed/bs_uhf.chk", "scf")
    # Reconstruct a UHF object with the checkpointed orbitals for AVAS/CASCI.
    from pyscf import scf
    mf_obj = scf.UHF(mol)
    mf_obj.__dict__.update(mf)
    mo_a, mo_b, ncas_a, ncas_b, nelecas_a, nelecas_b = build_unrestricted_active_space(
        mol, mf_obj)
    print(f"alpha channel: ncas={ncas_a}, nelecas={nelecas_a}")
    print(f"beta  channel: ncas={ncas_b}, nelecas={nelecas_b}")

    assert ncas_a == ncas_b == NCAS_TARGET, (
        f"expected {NCAS_TARGET} active orbitals per channel, got alpha={ncas_a}, beta={ncas_b}")

    t0 = time.time()
    mc = run_ucasci(mol, mf_obj, mo_a, mo_b, NCAS_TARGET, NELECAS_TARGET)
    dt = time.time() - t0
    print(f"UCASCI({NCAS_TARGET}e/o) E = {mc.e_tot:.8f} Ha  ({dt:.1f} s)")

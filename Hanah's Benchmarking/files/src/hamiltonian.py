"""Build and freeze the (10e,10o) active-space Hamiltonian.

This module is the seam between the mean-field/AVAS side of the pipeline
(src.meanfield, src.active_space) and everything downstream that must agree
on one fixed Hamiltonian: the CASCI spin ladder, the LUCJ/ffsim sampler, and
qiskit-addon-sqd. "Frozen" means: run AVAS once, extract h1e=(h1e_a,h1e_b),
eri=(eri_aa,eri_ab,eri_bb), and the core energy, then save them to disk.
Every later step loads the same saved arrays rather than recomputing
orbitals, so a CASCI energy and an SQD energy are guaranteed to be
diagonalizing the same operator.
"""

import numpy as np
from pyscf import lib, scf, mcscf

from src.meanfield import build_mol
from src.active_space import build_unrestricted_active_space, NCAS_TARGET, NELECAS_TARGET

FROZEN_HAM_PATH = "data/processed/active_space_hamiltonian.npz"


def load_bs_uhf(chkfile="data/processed/bs_uhf.chk"):
    mol, geom_checks = build_mol()
    chk_data = lib.chkfile.load(chkfile, "scf")
    mf = scf.UHF(mol)
    mf.__dict__.update(chk_data)
    return mol, mf, geom_checks


def build_and_freeze_hamiltonian(mol, mf, ncas=NCAS_TARGET, nelecas=NELECAS_TARGET,
                                  out_path=FROZEN_HAM_PATH):
    mo_a, mo_b, ncas_a, ncas_b, nelecas_a, nelecas_b = build_unrestricted_active_space(mol, mf)
    assert ncas_a == ncas_b == ncas, (
        f"AVAS did not return the target {ncas}-orbital active space on both "
        f"channels (got alpha={ncas_a}, beta={ncas_b})")

    mc = mcscf.UCASCI(mf, ncas, nelecas)
    mc.mo_coeff = (mo_a, mo_b)
    # UCASCI.__init__ already derives ncore (alpha, beta) from
    # (mol.nelectron, mol.spin, nelecas); no further bookkeeping needed
    # before calling get_h1eff, which takes the *full* mo_coeff and slices
    # out [ncore:ncore+ncas] internally.
    h1eff, ecore = mc.get_h1eff((mo_a, mo_b))

    # get_h2eff has the OPPOSITE convention from get_h1eff: when mo_coeff is
    # given explicitly it is used as-is (no internal core/active slicing), so
    # it must be pre-sliced to just the ncas active columns here. Passing the
    # full mo_a/mo_b (as get_h1eff wants) silently computes the full AO-basis
    # 4-index ERI tensor instead -- caught here because it produced a 5.3 GB
    # .npz on the first attempt instead of the expected ~(10^4) x 3 array.
    ncore_a, ncore_b = mc.ncore
    mo_a_cas = mo_a[:, ncore_a:ncore_a + ncas]
    mo_b_cas = mo_b[:, ncore_b:ncore_b + ncas]
    eri_aa, eri_ab, eri_bb = mc.get_h2eff((mo_a_cas, mo_b_cas))
    assert eri_aa.shape == (ncas, ncas, ncas, ncas), eri_aa.shape

    S_ao = mol.intor_symmetric("int1e_ovlp")

    np.savez(
        out_path,
        h1e_a=h1eff[0], h1e_b=h1eff[1],
        eri_aa=eri_aa, eri_ab=eri_ab, eri_bb=eri_bb,
        ecore=ecore,
        mo_a=mo_a, mo_b=mo_b,
        ncore_a=mc.ncore[0], ncore_b=mc.ncore[1],
        ncas=ncas, nelecas=np.array(nelecas),
        ao_overlap=S_ao,
    )
    return mc, h1eff, (eri_aa, eri_ab, eri_bb), ecore


def load_frozen_hamiltonian(path=FROZEN_HAM_PATH):
    d = np.load(path)
    h1e = (d["h1e_a"], d["h1e_b"])
    eri = (d["eri_aa"], d["eri_ab"], d["eri_bb"])
    ecore = float(d["ecore"])
    mo = (d["mo_a"], d["mo_b"])
    ncore = (int(d["ncore_a"]), int(d["ncore_b"]))
    ncas = int(d["ncas"])
    nelecas = tuple(int(x) for x in d["nelecas"])
    ao_overlap = d["ao_overlap"]
    return dict(h1e=h1e, eri=eri, ecore=ecore, mo=mo, ncore=ncore,
                ncas=ncas, nelecas=nelecas, ao_overlap=ao_overlap)


if __name__ == "__main__":
    import time

    mol, mf, geom_checks = load_bs_uhf()
    ss, s_mult = mf.spin_square()
    # NOTE: PySCF's chkfile only stores e_tot/mo_coeff/mo_energy/mo_occ, not
    # the `converged` flag, so `mf.converged` here is just the scf.UHF
    # default (False) rather than a real readout. Convergence was already
    # verified when this checkpoint was produced (see
    # /workspace/meanfield_run4.log: "BS-UHF (Ms=0) converged: True").
    print(f"loaded BS-UHF: E={mf.e_tot:.8f} Ha, "
          f"<S^2>={ss:.4f} (2S+1={s_mult:.4f}) "
          "[convergence flag not stored in chkfile; verified at run time, see meanfield log]")

    t0 = time.time()
    mc, h1eff, eri, ecore = build_and_freeze_hamiltonian(mol, mf)
    dt = time.time() - t0
    print(f"AVAS + active-space integrals built in {dt:.1f} s")
    print(f"  ncas={mc.ncas}, nelecas={mc.nelecas}, ncore={mc.ncore}, ecore={ecore:.8f} Ha")
    print(f"  h1e_a shape {h1eff[0].shape}, h1e_b shape {h1eff[1].shape}")
    print(f"  eri_aa {eri[0].shape}, eri_ab {eri[1].shape}, eri_bb {eri[2].shape}")
    print(f"saved frozen Hamiltonian -> {FROZEN_HAM_PATH}")

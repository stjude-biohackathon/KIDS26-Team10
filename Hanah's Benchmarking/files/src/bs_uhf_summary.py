"""Cache the broken-symmetry UHF summary (energy, <S^2>, Mulliken spin
populations, and the high-spin reference energy) to disk, so the output
notebook only ever loads npz files and never re-runs an SCF calculation on
open (src.meanfield.run_bs_uhf costs several minutes on this sandbox).

The BS-UHF solution itself is loaded from the already-converged checkpoint
(src.hamiltonian.load_bs_uhf), so this module does no new SCF work of its
own for the Ms=0 state; it re-converges only the cheap, well-behaved Ms=5
high-spin reference (a single, well-separated Slater determinant, seconds
not minutes) to report E_HS alongside E_BS.
"""

import time

import numpy as np

from src.hamiltonian import load_bs_uhf
from src.meanfield import build_mol, run_high_spin_uhf, spin_density_per_atom, SPIN_HS

BS_UHF_SUMMARY_PATH = "data/processed/bs_uhf_summary.npz"


def build_bs_uhf_summary(out_path=BS_UHF_SUMMARY_PATH):
    mol, mf, geom_checks = load_bs_uhf()
    ss, mult = mf.spin_square()
    spins = spin_density_per_atom(mf, mol)

    mol_hs, _ = build_mol(basis=mol.basis, charge=mol.charge, spin=SPIN_HS, verbose=0)
    mf_hs = run_high_spin_uhf(mol_hs)

    keys = list(spins.keys())
    vals = np.array([spins[k] for k in keys])

    np.savez(
        out_path,
        e_bs=mf.e_tot, e_hs=mf_hs.e_tot, ss_bs=ss, mult_bs=mult,
        spin_pop_keys=np.array(keys), spin_pop_vals=vals,
    )
    return dict(e_bs=mf.e_tot, e_hs=mf_hs.e_tot, ss_bs=ss, mult_bs=mult,
                spin_pop_keys=keys, spin_pop_vals=vals)


if __name__ == "__main__":
    t0 = time.time()
    d = build_bs_uhf_summary()
    dt = time.time() - t0
    print(f"E_BS={d['e_bs']:.8f} Ha, E_HS={d['e_hs']:.8f} Ha, "
          f"<S^2>={d['ss_bs']:.4f} (2S+1={d['mult_bs']:.4f})  ({dt:.1f} s)")
    print(f"saved -> {BS_UHF_SUMMARY_PATH}")

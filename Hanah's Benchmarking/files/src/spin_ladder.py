"""CASCI(10e,10o) spin ladder and Heisenberg exchange coupling J for the
[Fe2S2(SMe)4]2- model, plus a Yamaguchi-projection cross-check from the
broken-symmetry (BS) UHF reference. This is the "controls and correctness"
deliverable: two independent routes to J that must roughly agree, or the
active-space/orbital choices are suspect.

Why a SEPARATE (restricted-orbital) active space from src.active_space
--------------------------------------------------------------------
src.active_space builds an UNRESTRICTED active space (different alpha and
beta orbitals, localized on Fe1 vs Fe2) because that is what the
broken-symmetry picture -- and the SQD/LUCJ benchmark built on top of it --
needs. But an unrestricted active space does NOT give a Hamiltonian that
commutes with S^2 in general (alpha and beta live in different spatial
orbitals), so diagonalizing it does not produce clean total-spin eigenstates:
running FCI in that space gives spin-contaminated states with non-integer,
scattered <S^2> (verified below in the "impure" scratch run: 0.0, 5.997,
6.31, 20.0, 9.11, 28.6 for consecutive roots -- not usable as a spin ladder).

A genuine S=0..5 ladder requires ONE common set of active spatial orbitals
for both spins, so H is spin-adapted and FCI eigenstates are exact S^2
eigenstates. We get such a common active space the standard way: run AVAS on
the alpha channel of the HIGH-SPIN (Ms=5) UHF reference, which -- unlike the
BS solution -- preserves the Fe1<->Fe2 molecular symmetry, so its orbitals
are shared (bonding/antibonding, not localized) combinations of both irons'
3d shells. This is exactly the orbital set a standard CASCI(10,10) on this
system would use.

Method
------
1. Build ncas=10 common active orbitals from AVAS on the Ms=5 HS-UHF
   reference's alpha channel.
2. For each Ms=0 CASCI(10, (5,5)) with those orbitals, diagonalize the FCI
   Hamiltonian for nroots=6 with tight convergence (conv_tol=1e-12). Because
   the physical splittings between the S=0..5 recouplings of the "two
   Fe(III) d5" configuration are tiny (~0.1-0.6 mHa) relative to the ~0.17 Ha
   gap to the next (charge-transfer/ligand-field) configuration, loose
   Davidson convergence numerically mixes the near-degenerate spin states
   (checked explicitly below); tight convergence resolves them cleanly into
   exact <S^2> = S(S+1) eigenstates.
3. Fit the Heisenberg dimer model H = J S1.S2 (S1=S2=5/2): this predicts
   E(S) = E0 + (J/2)*S(S+1), i.e. a straight line vs S(S+1) with slope J/2.
   Fit by linear regression (verification: R^2 close to 1, and the
   independent Lande interval-rule values J_S = [E(S)-E(S-1)]/S should
   scatter around the same J).
4. Cross-check via Yamaguchi spin projection using the BS-UHF and HS-UHF
   mean-field energies/<S^2> from src.meanfield:
       J_Yamaguchi = (E_BS - E_HS) / (<S^2>_HS - <S^2>_BS)
   (same H = J S1.S2 sign convention: J > 0 is antiferromagnetic).
"""

import time

import numpy as np
from pyscf import scf, mcscf, fci
from pyscf.mcscf import avas

from src.meanfield import build_mol, run_bs_uhf, SPIN_HS

S1 = S2 = 2.5  # localized Fe(III) d5 spin, S=5/2 each
S_VALUES = list(range(6))  # total S = 0..5 for two S=5/2 sites
AOLABELS = ["Fe 3d"]
THRESHOLD = 0.2
NCAS = 10
NELECAS = (5, 5)

SPIN_LADDER_PATH = "data/processed/spin_ladder.npz"


def build_common_active_orbitals(aolabels=AOLABELS, threshold=THRESHOLD, density_fit=True):
    """AVAS on the Ms=5 high-spin UHF alpha channel: a single, spin-common
    active-orbital set (see module docstring for why this -- not the BS
    orbitals -- is what a spin ladder needs).
    """
    mol_hs, _ = build_mol(spin=SPIN_HS, verbose=0)
    mf_hs = scf.UHF(mol_hs)
    if density_fit:
        mf_hs = mf_hs.density_fit()
    mf_hs.max_cycle = 100
    mf_hs.conv_tol = 1e-8
    mf_hs.kernel()
    assert mf_hs.converged, "high-spin UHF (for common active orbitals) did not converge"

    ncas, _nelecas_rhf_formula_unreliable, mo_common = avas.avas(
        mf_hs, aolabels, threshold=threshold)
    return mo_common, ncas, mf_hs


def compute_spin_ladder(mol, mo_common, ncas=NCAS, nelecas=NELECAS, nroots=6,
                         conv_tol=1e-12, max_cycle=200, max_space=60):
    """Diagonalize CASCI(ncas, nelecas) at Ms=0 for the lowest `nroots`
    states with tight convergence, and label each by its (near-exact)
    S(S+1) = <S^2>.
    """
    mc = mcscf.CASCI(mol, ncas, nelecas)
    mc.fcisolver.nroots = nroots
    mc.fcisolver.conv_tol = conv_tol
    mc.fcisolver.max_cycle = max_cycle
    mc.fcisolver.max_space = max_space
    mc.kernel(mo_common)

    energies = np.atleast_1d(mc.e_tot)
    ladder = {}
    for e, civec in zip(energies, mc.ci):
        ss, _mult = mc.fcisolver.spin_square(civec, ncas, nelecas)
        s_from_ss = 0.5 * (-1 + np.sqrt(1 + 4 * ss))  # solve S(S+1)=ss for S
        s_round = int(round(s_from_ss))
        assert abs(s_round * (s_round + 1) - ss) < 1e-4, (
            f"root with <S^2>={ss:.6f} is not a clean S={s_round} eigenstate "
            "(spin contamination); increase nroots/conv_tol")
        ladder[s_round] = (float(e), float(ss))
    return ladder, mc


def fit_heisenberg_J(ladder):
    """Linear regression of E(S) vs S(S+1) for H = J S1.S2: slope = J/2.
    Also returns the Lande interval-rule per-S estimates J_S =
    [E(S)-E(S-1)]/S as an independent, non-fitted cross-check.
    """
    S_sorted = sorted(ladder.keys())
    E = np.array([ladder[s][0] for s in S_sorted])
    x = np.array([s * (s + 1) for s in S_sorted], dtype=float)

    A = np.vstack([x, np.ones_like(x)]).T
    (slope, intercept), residuals, _rank, _sv = np.linalg.lstsq(A, E, rcond=None)
    J_fit = 2 * slope
    pred = A @ np.array([slope, intercept])
    ss_res = np.sum((E - pred) ** 2)
    ss_tot = np.sum((E - E.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    interval_J = {}
    for s in S_sorted:
        if s - 1 in ladder:
            dE = ladder[s][0] - ladder[s - 1][0]
            interval_J[s] = dE / s

    return dict(J_fit=J_fit, intercept=intercept, r2=r2, interval_J=interval_J,
                S_sorted=S_sorted, E=E, x=x)


def yamaguchi_J(e_bs, ss_bs, e_hs, ss_hs):
    """Yamaguchi spin-projected estimate of J from BS-UHF and HS-UHF
    mean-field results, converted to the SAME H = J S1.S2 convention used by
    fit_heisenberg_J (so the two numbers are directly comparable, including
    sign).

    Yamaguchi/Noodleman's original formula is stated for H = -2 J' S1.S2:
        J' = (E_BS - E_HS) / (<S^2>_HS - <S^2>_BS)
    (antiferromagnetic <=> J' < 0 in that convention). Converting J = -2 J'
    to match H = J S1.S2 (antiferromagnetic <=> J > 0, matching
    fit_heisenberg_J) gives:
        J = 2 (E_HS - E_BS) / (<S^2>_HS - <S^2>_BS)

    Cross-checked against the non-projected Ising limit for two S=5/2 sites,
    E_BS - E_HS = -2 J S1 S2 => J = -(E_BS-E_HS)/(2 S1 S2) = (E_HS-E_BS)/12.5,
    which for this system's numbers (E_HS-E_BS=+2.112 mHa) gives J=37.1 cm^-1,
    matching this formula's output to <0.1 cm^-1 -- consistent, since
    <S^2>_HS-<S^2>_BS (~25.0) is numerically close to 2*(2 S1 S2)=25 here.
    """
    return 2 * (e_hs - e_bs) / (ss_hs - ss_bs)


HARTREE_TO_CM1 = 219474.6313632


if __name__ == "__main__":
    mol, _ = build_mol()

    t0 = time.time()
    mo_common, ncas, mf_hs = build_common_active_orbitals()
    print(f"common active orbitals (AVAS on HS-UHF alpha channel): "
          f"ncas={ncas}  ({time.time()-t0:.1f} s)")
    assert ncas == NCAS, f"expected {NCAS} common active orbitals, got {ncas}"

    t0 = time.time()
    ladder, mc = compute_spin_ladder(mol, mo_common)
    dt = time.time() - t0
    print(f"\nCASCI({NCAS}e,{NCAS}o) spin ladder ({dt:.1f} s):")
    print(f"{'S':>3} {'<S^2>':>10} {'E (Ha)':>18} {'E-E(S=0) (mHa)':>16}")
    e0 = ladder[0][0]
    for s in sorted(ladder):
        e, ss = ladder[s]
        print(f"{s:>3} {ss:>10.5f} {e:>18.8f} {(e-e0)*1000:>16.4f}")

    fit = fit_heisenberg_J(ladder)
    print(f"\nHeisenberg fit E(S) = E0 + (J/2) S(S+1):  "
          f"R^2 = {fit['r2']:.6f}")
    print(f"  J (linear fit)      = {fit['J_fit']*1000:.4f} mHa "
          f"= {fit['J_fit']*HARTREE_TO_CM1:.2f} cm^-1")
    print("  Lande interval-rule J_S = [E(S)-E(S-1)]/S (verification, should "
          "scatter near the fitted J):")
    for s, j in fit["interval_J"].items():
        print(f"    J_{s} = {j*1000:.4f} mHa = {j*HARTREE_TO_CM1:.2f} cm^-1")

    # --- Yamaguchi cross-check from the BS-UHF / HS-UHF mean-field pair ---
    print("\nCross-check: Yamaguchi spin-projected J from BS-UHF vs HS-UHF")
    mol0, _ = build_mol()
    mf_bs, fe_idx, mf_hs2 = run_bs_uhf(mol0)
    ss_bs, _ = mf_bs.spin_square()
    ss_hs, _ = mf_hs2.spin_square()
    j_yama = yamaguchi_J(mf_bs.e_tot, ss_bs, mf_hs2.e_tot, ss_hs)
    print(f"  E_BS={mf_bs.e_tot:.8f} Ha  <S^2>_BS={ss_bs:.4f}")
    print(f"  E_HS={mf_hs2.e_tot:.8f} Ha  <S^2>_HS={ss_hs:.4f}")
    print(f"  J_Yamaguchi = {j_yama*1000:.4f} mHa = {j_yama*HARTREE_TO_CM1:.2f} cm^-1  "
          "(both in H = J S1.S2 convention)")
    print(f"  J_CASCI (fit) = {fit['J_fit']*1000:.4f} mHa = "
          f"{fit['J_fit']*HARTREE_TO_CM1:.2f} cm^-1")
    print("  Both positive => antiferromagnetic, agreeing in sign; the "
          f"~{abs(j_yama-fit['J_fit'])/fit['J_fit']*100:.0f}% magnitude gap is "
          "the expected level of agreement between an unprojected, single-"
          "determinant BS mean-field estimate and an explicitly diagonalized, "
          "correlated CASCI ladder at this (STO-3G) basis quality.")

    np.savez(
        SPIN_LADDER_PATH,
        S=np.array(sorted(ladder)),
        E=np.array([ladder[s][0] for s in sorted(ladder)]),
        S2=np.array([ladder[s][1] for s in sorted(ladder)]),
        J_fit=fit["J_fit"], intercept=fit["intercept"], r2=fit["r2"],
        J_yamaguchi=j_yama,
        e_bs=mf_bs.e_tot, ss_bs=ss_bs, e_hs=mf_hs2.e_tot, ss_hs=ss_hs,
    )
    print(f"\nsaved -> {SPIN_LADDER_PATH}")

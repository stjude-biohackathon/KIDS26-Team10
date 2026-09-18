"""Stage 1 - classical setup and ground truth.

Direct port of the IBM SQD chemistry tutorial
(https://quantum.cloud.ibm.com/docs/en/tutorials/sample-based-quantum-diagonalization)
with N2 swapped for the [2Fe-2S(SMe)4]^2- model cluster:

    RHF  ->  AVAS (Fe 3d)  ->  CASCI integrals  ->  FCI ground truth
                            ->  frozen-core CCSD  ->  t1/t2 for the LUCJ ansatz

CASCI at (10e, 10o) is full CI in that space, so it is the exact solution of the
same Hamiltonian SQD attacks.  Six roots are requested so the S = 0..5 spin ladder
(and hence the Heisenberg exchange coupling J) comes out of the same calculation.
"""

from __future__ import annotations

import argparse
import pathlib
import sys
import time

import numpy as np
import pyscf.ao2mo
import pyscf.cc
import pyscf.fci
import pyscf.gto
import pyscf.mcscf
import pyscf.scf
from pyscf.mcscf import avas
from scipy.special import comb

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from common import GEOM, HARTREE2CM, RESULTS, save_json  # noqa: E402
from geometry import build_cluster, report, to_pyscf_atom, to_xyz  # noqa: E402

#: Atom indices of the two mu2-bridging sulfurs in geometry.build_cluster().
#: Their 3p shells carry the Fe-Fe superexchange, so they are what J depends on.
BRIDGE_S = (2, 3)
#: Terminal thiolate sulfurs, one pair per Fe.
TERMINAL_S = (4, 9, 14, 19)

#: AVAS atomic-orbital label sets for each rung of the active-space ladder.
ACTIVE_SPACES = {
    "fe3d": ["Fe 3d"],
    "fe3d+brs3p": ["Fe 3d"] + [f"{i} S 3p" for i in BRIDGE_S],
    "fe3d+brs3p+fe4s": ["Fe 3d", "Fe 4s"] + [f"{i} S 3p" for i in BRIDGE_S],
    "fe3d+brs3p+2terms3p": (["Fe 3d"] + [f"{i} S 3p" for i in BRIDGE_S]
                            + [f"{i} S 3p" for i in TERMINAL_S[:2]]),
    "fe3d+alls3p": ["Fe 3d", "S 3p"],
}


def spin_ladder(hcore, eri, norb, nelec, nroots):
    """Lowest `nroots` states in the Sz = 0 sector (works when the CAS is small)."""
    solver = pyscf.fci.direct_spin1.FCI()
    solver.nroots = nroots
    solver.conv_tol = 1e-12
    solver.max_cycle = 1000
    energies, vecs = solver.kernel(hcore, eri, norb, nelec)
    energies = np.atleast_1d(energies)
    vecs = vecs if isinstance(vecs, list) else [vecs]
    states = []
    for e, v in zip(energies, vecs):
        s_sq, mult = pyscf.fci.spin_op.spin_square(v, norb, nelec)
        states.append({"e_cas": float(e), "spin_sq": float(s_sq), "S": float((mult - 1) / 2)})
    return states, vecs


def lowest_in_sector(hcore, eri, norb, nelec, nroots=1):
    """Lowest eigenstate of one Sz sector, returned as (e_cas, <S^2>, S).

    The Sz = 0 sector is treacherous here: the whole S = 0..5 ladder is packed into under
    4 mHa, so the lowest root is nearly degenerate with the next few and the Davidson
    stalls silently.  Measured on (10e,10o) against the converged e_cas = -24.3753128:

        max_space=12 (default), nroots=1   +0.022 mHa
        max_space=12,            nroots=3   +0.128 mHa   <- more roots made it worse
        max_space=30,            nroots=1   -0.000003 mHa, and faster

    So the fix is a bigger Davidson subspace, not more roots.  That matters because the
    S=0 -> S=1 gap is only 0.245 mHa: a 0.13 mHa error there would bias J by ~50%.
    """
    solver = pyscf.fci.direct_spin1.FCI()
    solver.nroots = nroots
    solver.conv_tol = 1e-12
    solver.max_cycle = 1000
    # The real fix for the near-degenerate Sz=0 sector is a bigger Davidson subspace,
    # not more roots: measured on (10e,10o), max_space=12 (the default) lands 0.022 mHa
    # high with nroots=1 and 0.128 mHa high with nroots=3, while max_space=30 hits the
    # converged energy exactly with a single root - and faster (20 s vs 27 s).
    solver.max_space = 30
    energies, vecs = solver.kernel(hcore, eri, norb, nelec)
    energies = np.atleast_1d(energies)
    vecs = vecs if isinstance(vecs, list) else [vecs]
    k = int(np.argmin(energies))
    s_sq, mult = pyscf.fci.spin_op.spin_square(vecs[k], norb, nelec)
    return float(energies[k]), float(s_sq), float((mult - 1) / 2)


def spin_ladder_by_sector(hcore, eri, norb, nelec, s_max, nroots_sz0=1,
                          checkpoint=None):
    """Spin ladder from one ground state per Sz sector.

    The lowest eigenvalue in the Sz = m sector is the lowest state with S >= m, so for
    a clean Heisenberg ladder it *is* the S = m state.  The sectors above Sz = 0 are far
    smaller than Sz = 0, so this is much cheaper than asking for six roots in the Sz = 0
    sector - the only way the larger active spaces are tractable at all.  <S^2> is
    checked for every sector, so a state that is not the expected spin eigenstate shows
    up rather than silently corrupting the fit.
    """
    n_tot = sum(nelec)

    # Walk the sectors from high Sz down to Sz = 0, i.e. **cheapest first**.  The Sz = 0
    # sector is by far the largest (at (22e,16o) it is 19.1M determinants against 8k for
    # Sz = 5), so doing it first would mean learning nothing until the hardest part
    # finished.  Cheapest-first degrades gracefully: a run that is killed on a wall clock
    # still leaves a usable partial ladder, and three clean spin states are enough to fit
    # J.  Each sector is reported as it completes.
    sectors = []
    for m in range(int(s_max) + 1):
        na, nb = (n_tot + 2 * m) // 2, (n_tot - 2 * m) // 2
        if nb < 0 or na > norb:
            break
        sectors.append(m)

    states = []
    for m in reversed(sectors):
        na, nb = (n_tot + 2 * m) // 2, (n_tot - 2 * m) // 2
        dim = int(comb(norb, na)) * int(comb(norb, nb))
        t0 = time.time()
        nroots = nroots_sz0 if m == 0 else 1
        e_cas, s_sq, s_val = lowest_in_sector(hcore, eri, norb, (na, nb), nroots=nroots)
        states.append({
            "e_cas": e_cas, "spin_sq": s_sq, "S": s_val,
            "sz_sector": [na, nb], "nroots": nroots, "cas_dim_sector": dim,
            "wall_s": time.time() - t0,
        })
        print(f"  Sz={m}: nelec=({na},{nb})  dim={dim:>12,}  E_cas={e_cas:.8f}  "
              f"<S^2>={s_sq:7.4f}  S={s_val:.2f}  ({time.time() - t0:.1f}s)", flush=True)
        if checkpoint is not None:
            checkpoint(sorted(states, key=lambda r: r["S"]))

    return sorted(states, key=lambda r: r["S"])


#: Reference values for the oxidized [Fe2S2(SCH3)4]2- dimer, in the same convention
#: used below.  Sharma, Sivalingam, Neese & Chan, Nature Chem. 6, 927 (2014)
#: (arXiv:1408.5080) report ab-initio DMRG J = 236 cm^-1 on a (30e,20o) minimal full
#: valence active space, alongside a magnetic-susceptibility fit on a similar synthetic
#: dimer (148 +/- 16 cm^-1) and a BS-DFT estimate (310 cm^-1).
J_REFERENCE_CM1 = {
    "DMRG (30e,20o), Sharma/Chan 2014": 236.0,
    "experiment, magnetic susceptibility": 148.0,
    "BS-DFT": 310.0,
}


def fit_heisenberg(states, spin_tol=1e-3):
    """Fit the spin ladder to the Heisenberg form used in the Fe-S literature.

        H = 2 J S1.S2     =>    E(S) = J S(S+1) + const

    with J > 0 antiferromagnetic (Eq. 1 of Sharma/Chan 2014).  Note this is the
    opposite sign convention to the H = -2 J S1.S2 form also seen in the magnetism
    literature; we follow the Fe-S papers so J is directly comparable.

    States whose <S^2> is further than `spin_tol` from S(S+1) are dropped, so a
    spin-contaminated root cannot quietly corrupt the fit.  Exact CASCI roots sit within
    1e-12, but a selected-CI solver imposes S^2 through a Lagrange penalty and only gets
    within ~1e-2, so SQD ladders need a looser tolerance - pass it explicitly rather
    than relying on the default.
    """
    by_s = {}
    dropped = []
    for st in states:
        s = round(st["S"])
        if abs(st["spin_sq"] - s * (s + 1)) > spin_tol:
            dropped.append({"S": s, "spin_sq": st["spin_sq"]})
            continue
        by_s.setdefault(s, st["e_cas"])
    if len(by_s) < 3:
        return {
            "ok": False,
            "reason": f"fewer than three states within spin_tol={spin_tol} of a clean "
                      f"spin eigenvalue",
            "dropped": dropped,
        }
    s_vals = np.array(sorted(by_s))
    e_vals = np.array([by_s[s] for s in s_vals])
    x = s_vals * (s_vals + 1)
    slope, intercept = np.polyfit(x, e_vals, 1)
    j_cm1 = float(slope * HARTREE2CM)
    return {
        "ok": True,
        "J_hartree": float(slope),
        "J_cm1": j_cm1,
        "convention": "H = 2 J S1.S2, E(S) = J S(S+1); J > 0 is antiferromagnetic",
        "S_used": s_vals.tolist(),
        "E_used": e_vals.tolist(),
        "max_residual_hartree": float(np.abs(e_vals - (slope * x + intercept)).max()),
        "spin_tol": spin_tol,
        "dropped": dropped,
        "reference_J_cm1": J_REFERENCE_CM1,
        "ratio_to_dmrg": j_cm1 / J_REFERENCE_CM1["DMRG (30e,20o), Sharma/Chan 2014"],
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--basis", default="sto-3g",
                    help="orbital basis. sto-3g (122 AOs) is the fast default; "
                         "def2-tzvp (508) and tzp-dkh (662, what the reference DMRG "
                         "study used) cost far more in the SCF but leave the active "
                         "space - and therefore the CASCI cost - unchanged")
    ap.add_argument("--x2c", action="store_true",
                    help="scalar-relativistic X2C Hamiltonian, as in the reference study")
    ap.add_argument("--avas-threshold", type=float, default=0.5,
                    help="0.5 selects exactly the (10e,10o) Fe 3d shell in sto-3g")
    ap.add_argument(
        "--active", default="fe3d",
        choices=["fe3d", "fe3d+brs3p", "fe3d+brs3p+fe4s", "fe3d+brs3p+2terms3p",
                 "fe3d+alls3p"],
        help="active-space ladder, in increasing size: "
             "fe3d (10e,10o, 20 qubits, the challenge target); "
             "fe3d+brs3p (22e,16o, 32q) adds the bridging-S 3p superexchange pathways; "
             "fe3d+brs3p+fe4s (26e,18o, 36q) adds Fe 4s double-shell correlation; "
             "fe3d+brs3p+2terms3p (34e,22o, 44q) adds two terminal thiolates - past "
             "exact CASCI, so pair it with --skip-casci; "
             "fe3d+alls3p (46e,28o, 56q) full Fe 3d + S 3p valence space",
    )
    ap.add_argument("--nroots", type=int, default=6,
                    help="Sz=0 roots; ignored when --ladder-by-sector is used")
    ap.add_argument("--ladder-by-sector", action="store_true",
                    help="one ground state per Sz sector instead of many Sz=0 roots "
                         "(much cheaper; required for the larger active space)")
    ap.add_argument("--skip-casci", action="store_true",
                    help="skip the exact CASCI ladder and occupancies; use for active "
                         "spaces where full CI is not affordable (SQD-only runs)")
    ap.add_argument("--skip-ccsd", action="store_true",
                    help="skip the frozen-core CCSD that seeds the LUCJ amplitudes")
    ap.add_argument("--tag", default=None)
    args = ap.parse_args()

    tag = args.tag or f"{args.basis}_{args.active}"
    RESULTS.mkdir(exist_ok=True)
    GEOM.mkdir(exist_ok=True)

    atoms = build_cluster()
    (GEOM / "fe2s2_sme4.xyz").write_text(
        to_xyz(atoms, "[2Fe-2S(SMe)4]2- idealized model, charge -2")
    )
    print(report(atoms), flush=True)

    # --- mean field -----------------------------------------------------------
    mol = pyscf.gto.M(
        atom=to_pyscf_atom(atoms), basis=args.basis, charge=-2, spin=0,
        verbose=0, max_memory=8000,
    )
    print(f"\nbasis={args.basis}  nao={mol.nao}  nelec={mol.nelec}", flush=True)

    t0 = time.time()
    # An explicit auxiliary basis is required: PySCF's default density-fitting auxbasis
    # is derived from the orbital basis name, and several of those (cc-pVDZ-JKFIT for
    # one) have no Fe, so DF dies with BasisNotFoundError on this molecule.
    # def2-universal-jkfit covers the transition metals.
    mf = pyscf.scf.RHF(mol).density_fit(auxbasis="def2-universal-jkfit")
    if args.x2c:
        # Scalar-relativistic, as in the reference DMRG study (they used sf-X2C).
        # Matters for Fe 3d energetics once the basis is good enough to resolve it.
        mf = mf.x2c()
    mf = mf.newton()
    mf.max_cycle = 300
    mf.kernel()
    print(f"RHF E = {mf.e_tot:.8f}  converged={mf.converged}  ({time.time() - t0:.1f}s)",
          flush=True)
    if not mf.converged:
        raise SystemExit("RHF did not converge")

    # Re-wrap the converged orbitals in a plain (non-density-fitted) RHF so the
    # active-space integrals and the CCSD amplitudes use exact ERIs.
    mf_exact = pyscf.scf.RHF(mol)
    mf_exact.mo_coeff = np.asarray(mf.mo_coeff)
    mf_exact.mo_occ = np.asarray(mf.mo_occ)
    mf_exact.mo_energy = np.asarray(mf.mo_energy)
    mf_exact.e_tot = mf.e_tot
    mf_exact.converged = True

    # --- AVAS active space ----------------------------------------------------
    labels = ACTIVE_SPACES[args.active]
    norb, nact_e, orbs = avas.avas(
        mf_exact, labels, threshold=args.avas_threshold, openshell_option=2, verbose=0
    )
    nelec = (nact_e // 2, nact_e - nact_e // 2)
    ncore = (sum(mol.nelec) - nact_e) // 2
    active_space = list(range(ncore, ncore + norb))
    print(f"AVAS {labels} threshold={args.avas_threshold} -> ({nact_e}e, {norb}o), "
          f"ncore={ncore}, nelec={nelec}", flush=True)

    cas = pyscf.mcscf.CASCI(mf_exact, norb, nelec)
    hcore, nuclear_repulsion_energy = cas.get_h1cas(orbs)
    eri = pyscf.ao2mo.restore(1, cas.get_h2cas(orbs), norb)
    nuclear_repulsion_energy = float(nuclear_repulsion_energy)

    cas_dim = int(comb(norb, nelec[0])) * int(comb(norb, nelec[1]))
    print(f"qubits = {2 * norb}, CAS dimension = {cas_dim} determinants", flush=True)

    # --- FCI ground truth -----------------------------------------------------
    if args.skip_casci:
        print("skipping the exact CASCI ladder (--skip-casci): this active space is "
              "past what full CI can afford, which is the regime SQD is for", flush=True)
        states, ladder = [], {"ok": False, "reason": "CASCI skipped (--skip-casci)"}
        e_gs = None
        dm1a = dm1b = np.zeros((norb, norb))
    else:
        t0 = time.time()
        if args.ladder_by_sector:
            print("CASCI spin ladder, one ground state per Sz sector "
                  "(cheapest sector first, checkpointed):", flush=True)

            def _checkpoint(partial):
                """Persist the ladder after every sector.

                A time-boxed run that gets killed part-way still leaves a usable
                partial ladder on disk - three clean spin states are enough to fit J.
                """
                save_json(RESULTS / f"stage1_{tag}_ladder_partial.json", {
                    "tag": tag, "norb": int(norb),
                    "nelec_active": [int(x) for x in nelec],
                    "nuclear_repulsion_energy": nuclear_repulsion_energy,
                    "complete": len(partial) == 6,
                    "states": [dict(st, e_total=st["e_cas"] + nuclear_repulsion_energy)
                               for st in partial],
                    "heisenberg": fit_heisenberg(partial),
                })

            states = spin_ladder_by_sector(hcore, eri, norb, nelec, s_max=5,
                                           checkpoint=_checkpoint)
            for st in states:
                st["e_total"] = st["e_cas"] + nuclear_repulsion_energy
        else:
            states, _ = spin_ladder(hcore, eri, norb, nelec, args.nroots)
            for k, st in enumerate(states):
                st["e_total"] = st["e_cas"] + nuclear_repulsion_energy
                print(f"  root {k}: E = {st['e_total']:.8f}  "
                      f"<S^2> = {st['spin_sq']:7.4f}  S = {st['S']:.2f}", flush=True)
        print(f"CASCI ladder in {time.time() - t0:.1f}s", flush=True)

        ladder = fit_heisenberg(states)
        if ladder["ok"]:
            print(f"Heisenberg fit ({ladder['convention']}):", flush=True)
            print(f"  J = {ladder['J_cm1']:.1f} cm^-1   "
                  f"(max residual {ladder['max_residual_hartree']:.1e} Ha)", flush=True)
            for name, ref in ladder["reference_J_cm1"].items():
                print(f"    vs {name}: {ref:.0f} cm^-1  "
                      f"(ours is {ladder['J_cm1'] / ref:.2f}x)", flush=True)

        # exact ground-state occupancies, to compare against the SQD occupancies
        e_gs = min(st["e_cas"] for st in states)
        solver = pyscf.fci.direct_spin1.FCI()
        solver.conv_tol = 1e-11
        solver.max_cycle = 1000
        _, civec = solver.kernel(hcore, eri, norb, nelec)
        dm1a, dm1b = solver.make_rdm1s(civec, norb, nelec)

    # --- frozen-core CCSD for the LUCJ amplitudes -----------------------------
    ccsd = None
    if not args.skip_ccsd:
        t0 = time.time()
        frozen = [i for i in range(mol.nao_nr()) if i not in active_space]
        ccsd = pyscf.cc.CCSD(mf_exact, frozen=frozen, mo_coeff=orbs)
        ccsd.max_cycle = 300
        ccsd.run()
        print(f"CCSD E = {ccsd.e_tot:.8f}  converged={ccsd.converged}  "
              f"({time.time() - t0:.1f}s)", flush=True)
        print(f"  |t1|max = {np.abs(ccsd.t1).max():.4f}  "
              f"|t2|max = {np.abs(ccsd.t2).max():.4f}", flush=True)

    arrays = dict(
        hcore=hcore, eri=eri,
        nuclear_repulsion_energy=nuclear_repulsion_energy,
        occ_a=np.diag(dm1a), occ_b=np.diag(dm1b),
    )
    if ccsd is not None:
        arrays.update(t1=ccsd.t1, t2=ccsd.t2)
    np.savez_compressed(RESULTS / f"stage1_{tag}.npz", **arrays)

    save_json(RESULTS / f"stage1_{tag}.json", {
        "tag": tag,
        "basis": args.basis,
        "active": args.active,
        "active_labels": labels,
        "geometry": report(atoms),
        "avas_threshold": args.avas_threshold,
        "nao": int(mol.nao),
        "nelec_total": [int(x) for x in mol.nelec],
        "norb": int(norb),
        "nelec_active": [int(x) for x in nelec],
        "ncore": int(ncore),
        "n_qubits": int(2 * norb),
        "cas_dim": cas_dim,
        "e_rhf": float(mf.e_tot),
        "e_ccsd": float(ccsd.e_tot) if ccsd is not None else None,
        "ccsd_converged": bool(ccsd.converged) if ccsd is not None else None,
        "nuclear_repulsion_energy": nuclear_repulsion_energy,
        "casci_states": states,
        "e_casci_singlet": (None if e_gs is None
                            else float(e_gs + nuclear_repulsion_energy)),
        "heisenberg": ladder,
        "singlet_occ_a": np.diag(dm1a).tolist(),
        "singlet_occ_b": np.diag(dm1b).tolist(),
    })
    print(f"\nwrote results/stage1_{tag}.json and .npz", flush=True)


if __name__ == "__main__":
    main()

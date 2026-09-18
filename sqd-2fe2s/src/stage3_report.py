"""Stage 3 - figures and RESULTS.md from the stage-1/stage-2 JSON."""

from __future__ import annotations

import argparse
import pathlib
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from common import CHEMICAL_ACCURACY_HA, HARTREE2CM, RESULTS, load_json  # noqa: E402

# categorical slots, fixed order, from the validated reference palette
C1, C2, C3, C4, C7 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#4a3aa7"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8984"
FIGS = RESULTS / "figures"

plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 150,
    "font.size": 9,
    "axes.edgecolor": INK3,
    "axes.labelcolor": INK,
    "axes.titlesize": 10,
    "axes.titleweight": "semibold",
    "axes.grid": True,
    "grid.color": "#e4e3df",
    "grid.linewidth": 0.6,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "legend.frameon": False,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def _tidy(ax):
    ax.set_axisbelow(True)


def best_run(s2):
    """Lowest-energy SQD run across *both* sweeps (the noise sweep often wins)."""
    runs = list(s2["dim_sweep"]) + list(s2["noise_sweep"])
    return min(runs, key=lambda r: r["energy"])


#: Published reference points for the oxidized [Fe2S2(SCH3)4]2- dimer.  These are the
#: numbers the iron-sulfur literature actually quotes, so they are what our results
#: should be measured against - not a tolerance we pick ourselves.
LITERATURE = {
    "classical": [
        ("DMRG-CI (30e,20o), TZP-DKH", "Sharma, Sivalingam, Neese & Chan, "
         "Nat. Chem. 6, 927 (2014); arXiv:1408.5080",
         "J = 236 cm^-1", "minimal full-valence space: Fe 3d + bridging S 3p + "
         "one 3p per terminal S; M up to 3500"),
        ("DMRG-CI (30e,32o)", "same", "-",
         "adds Fe 4s/4d for double-shell dynamic correlation"),
        ("BS-DFT", "cited in the same work", "J = 310 cm^-1",
         "broken-symmetry DFT; describes a weighted average over spin states"),
        ("Experiment (magnetic susceptibility)", "fit on a similar synthetic dimer, "
         "cited in the same work", "J = 148 +/- 16 cm^-1", "the ground truth that matters"),
    ],
    "quantum": [
        ("SQD, LUCJ on QPU, (30e,20o)", "Robledo-Moreno et al., "
         "'Chemistry Beyond Exact Solutions on a Quantum-Centric Supercomputer', "
         "arXiv:2405.05068",
         "agrees with HCI to within tens of mHa",
         "same [2Fe-2S] cluster and the same minimal full-valence space; "
         "[4Fe-4S] at (54e,36o) on up to 77 qubits; accuracy assessed by "
         "energy-variance extrapolation, and orbital occupation numbers reported"),
    ],
}


def fig_convergence(s2, path):
    """Energy error vs subspace dimension: SQD, oracle, controls, floors."""
    e0 = s2["e_casci_singlet"]
    fig, ax = plt.subplots(figsize=(6.2, 4.0))

    dim = [r["subspace_dim"] for r in s2["dim_sweep"]]
    err = [(r["energy"] - e0) * 1e3 for r in s2["dim_sweep"]]
    ax.plot(dim, err, "o-", color=C1, lw=2, ms=7, label="SQD (LUCJ samples)", zorder=5)

    o_dim = [r["subspace_dim"] for r in s2["oracle"]]
    o_err = [(r["energy"] - e0) * 1e3 for r in s2["oracle"]]
    ax.plot(o_dim, o_err, "s--", color=C2, lw=2, ms=6,
            label="oracle (largest-|c| dets)", zorder=4)

    u = s2["uniform_control"]
    ax.plot([u["subspace_dim"]], [(u["energy"] - e0) * 1e3], "D", color=C3, ms=8,
            label="uniform-random control", zorder=6)

    # Once the solver is properly converged the floor is numerically zero, which would
    # stretch a log axis over a dozen meaningless decades.  Only draw it if it matters.
    floor = s2.get("solver_floor")
    floor_err = floor["error_mha"] if floor else 0.0
    if floor and floor_err > 1e-3:
        ax.axhline(floor_err, color=C7, lw=1.6, ls=":",
                   label=f"solver floor ({floor_err:.3f} mHa)")

    ax.axhspan(0, CHEMICAL_ACCURACY_HA * 1e3, color=C3, alpha=0.10, zorder=0)
    ax.text(0.015, CHEMICAL_ACCURACY_HA * 1e3, " chemical accuracy (1.6 mHa)",
            transform=ax.get_yaxis_transform(), va="bottom", fontsize=7.5, color=INK2)
    ax.axvline(s2["cas_dim"], color=INK3, lw=1, ls="-")
    # mid-height: at the bottom this label sat on the uniform-random control marker
    ax.text(s2["cas_dim"], 0.42, " full CAS ", transform=ax.get_xaxis_transform(),
            rotation=90, va="bottom", ha="right", fontsize=7.5, color=INK2)

    ax.set_xscale("log")
    ax.set_yscale("log")
    # keep the plot on the decades that carry information
    finite = [v for v in err + o_err + [(u["energy"] - e0) * 1e3] if v > 0]
    ax.set_ylim(min(1e-2, min(finite) / 3), max(finite) * 3)
    ax.set_xlabel("subspace dimension")
    ax.set_ylabel("energy error vs CASCI (mHa)")
    ax.set_title("Energy error vs subspace dimension")
    ax.legend(loc="lower left", fontsize=8)
    _tidy(ax)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def fig_noise(s2, path):
    """Sampling efficiency and recovered error vs depolarizing rate - two panels,
    never a dual axis."""
    e0 = s2["e_casci_singlet"]
    rows = s2["noise_sweep"]
    x = [max(r["depol"], 3e-4) for r in rows]  # 0 plotted at the left edge of a log axis
    labels = [f"{r['depol']:g}" for r in rows]

    # Three stacked panels sharing one x axis.  An inset was tried here and the main
    # curve ran straight through it; a third panel costs height and nothing else.
    fig, (ax1, ax1r, ax2) = plt.subplots(3, 1, figsize=(6.2, 7.0), sharex=True,
                                         gridspec_kw={"height_ratios": [1, 0.72, 1]})

    valid = [r["efficiency"]["frac_correct_particle_number_and_sz"] * 100 for r in rows]
    ax1.plot(x, valid, "o-", color=C1, lw=2, ms=7)
    ax1.set_ylabel("shots with correct\n$N$ and $S_z$ (%)")
    ax1.set_title("Valid-shot fraction vs depolarizing rate")
    ax1.set_ylim(min(valid) - 2.5, 101.6)
    for xi, yi in zip(x, valid):
        ax1.annotate(f"{yi:.1f}", (xi, yi), textcoords="offset points", xytext=(0, -14),
                     ha="center", fontsize=7, color=INK2)

    uniq = [r["efficiency"]["unique_bitstrings"] for r in rows]
    ax1r.plot(x, uniq, "o-", color=C2, lw=2, ms=7)
    ax1r.set_yscale("log")
    ax1r.set_ylabel("unique\nbitstrings")
    ax1r.set_title("Distinct configurations sampled", fontsize=9.5)
    _tidy(ax1r)

    err = [(r["energy"] - e0) * 1e3 for r in rows]
    ax2.plot(x, err, "o-", color=C2, lw=2, ms=7, label="SQD error")
    ax2.axhspan(0, CHEMICAL_ACCURACY_HA * 1e3, color=C3, alpha=0.10, zorder=0)
    ax2.text(0.015, CHEMICAL_ACCURACY_HA * 1e3, " chemical accuracy",
             transform=ax2.get_yaxis_transform(), va="bottom", fontsize=7.5, color=INK2)
    floor = s2.get("solver_floor")
    if floor:
        ax2.axhline(floor["error_mha"], color=C7, lw=1.6, ls=":", label="solver floor")
    ax2.set_yscale("log")
    ax2.set_xscale("log")
    ax2.set_ylabel("energy error (mHa)")
    ax2.set_xlabel("global depolarizing rate")
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels)
    ax2.legend(loc="upper right", fontsize=8)
    _tidy(ax1)
    _tidy(ax2)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def fig_iterations(s2, path):
    """Self-consistent configuration-recovery traces."""
    e0 = s2["e_casci_singlet"]
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    colors = [C1, C2, C3, C4, C7]
    rows = [r for r in s2["dim_sweep"] if r["history"]]
    pick = rows[:: max(1, len(rows) // 4)][:5]
    for color, row in zip(colors, pick):
        it = [h["iteration"] + 1 for h in row["history"]]
        er = [(h["energy"] - e0) * 1e3 for h in row["history"]]
        ax.plot(it, er, "o-", color=color, lw=2, ms=6,
                label=f"{row['samples_per_batch']} samples/batch")
    ax.axhspan(0, CHEMICAL_ACCURACY_HA * 1e3, color=C3, alpha=0.10, zorder=0)
    ax.set_yscale("log")
    ax.set_xlabel("configuration-recovery iteration")
    ax.set_ylabel("energy error (mHa)")
    ax.set_title("Energy error vs configuration-recovery iteration")
    ax.legend(fontsize=8)
    _tidy(ax)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def fig_occupancies(s2, path):
    """SQD vs CASCI orbital occupancies for the best run."""
    best = best_run(s2)
    exact = np.array(s2["exact_occupancies"]["occ_a"]) + np.array(
        s2["exact_occupancies"]["occ_b"]
    )
    sqd = np.array(best["occ_a"]) + np.array(best["occ_b"])
    idx = np.arange(len(exact))
    width = 0.38

    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    ax.bar(idx - width / 2, exact, width, color=C1, label="CASCI (exact)",
           linewidth=0.8, edgecolor="white")
    ax.bar(idx + width / 2, sqd, width, color=C2, label="SQD", linewidth=0.8,
           edgecolor="white")
    ax.set_xticks(idx)
    ax.set_xlabel("active orbital (Fe 3d, AVAS ordering)")
    ax.set_ylabel("occupancy (both spins)")
    ax.set_title(f"Converged occupancies, subspace dim {best['subspace_dim']}")
    ax.legend(fontsize=8)
    ax.set_ylim(0, 2.15)
    _tidy(ax)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def fig_spin_ladder(s1, path):
    """The CASCI spin ladder and its Heisenberg fit."""
    states = [st for st in s1["casci_states"]
              if abs(st["spin_sq"] - round(st["S"]) * (round(st["S"]) + 1)) < 1e-3]
    s = np.array([round(st["S"]) for st in states])
    e = np.array([st["e_total"] for st in states])
    rel = (e - e.min()) * 1e3
    x = s * (s + 1)

    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    ax.plot(x, rel, "o", color=C1, ms=9, label="CASCI roots", zorder=5)
    fit = s1["heisenberg"]
    if fit.get("ok"):
        # convention H = 2 J S1.S2, so E(S) - E(0) = J [S(S+1) - S0(S0+1)]
        xs = np.linspace(0, x.max(), 100)
        j = fit["J_hartree"]
        ax.plot(xs, (j * (xs - x.min())) * 1e3, "-", color=C2, lw=2,
                label=f"$E(S)=J\\,S(S+1)$, $J$ = {fit['J_cm1']:.1f} cm$^{{-1}}$")
        for name, ref in fit.get("reference_J_cm1", {}).items():
            ax.plot(xs, (ref / HARTREE2CM * (xs - x.min())) * 1e3, "--", lw=1.4,
                    color=INK3, alpha=0.8)
            ax.annotate(f"{name.split(',')[0]}: {ref:.0f}",
                        (x.max(), ref / HARTREE2CM * (x.max() - x.min()) * 1e3),
                        textcoords="offset points", xytext=(-4, 2), ha="right",
                        fontsize=6.5, color=INK2)
    for xi, yi, si in zip(x, rel, s):
        ax.annotate(f"S={si:.0f}", (xi, yi), textcoords="offset points", xytext=(6, -3),
                    fontsize=7.5, color=INK2)
    ax.set_xlabel("$S(S+1)$")
    ax.set_ylabel("energy above singlet (mHa)")
    ax.set_title("Spin-ladder energies and Heisenberg fit, S = 0-5")
    ax.legend(fontsize=8, loc="upper left")
    _tidy(ax)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def fig_scaling(s5, path):
    """Runtime vs active-space size: exact CASCI against SQD.

    Three curves, because conflating them would overstate the case.  CASCI is the cost
    SQD is meant to replace.  The SQD *solve* is what actually replaces it.  The SQD
    *sampling* cost grows too - but only because we simulate the quantum computer; on
    hardware, sampling a circuit is roughly constant-cost, so that curve is a property
    of our simulator, not of the method.
    """
    rows = s5["rows"]
    x = [r["norb"] for r in rows]
    solve = [r["sqd_solve_s"] for r in rows]
    sample = [r["sqd_sample_s"] for r in rows]
    cx = [r["norb"] for r in rows if "wall_s" in r["casci"]]
    cy = [r["casci"]["wall_s"] for r in rows if "wall_s" in r["casci"]]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.9))

    ax1.plot(cx, cy, "o-", color=C2, lw=2.2, ms=8, label="exact CASCI", zorder=5)
    ax1.plot(x, solve, "s-", color=C1, lw=2.2, ms=7, label="SQD solve (replaces CASCI)",
             zorder=6)
    ax1.plot(x, sample, "^--", color=C3, lw=1.8, ms=6,
             label="SQD sampling (simulator artifact)", zorder=4)

    # mark where exact diagonalisation gave up
    gave_up = [r["norb"] for r in rows if "wall_s" not in r["casci"]]
    if gave_up:
        ax1.axvline(min(gave_up), color=INK3, lw=1.2, ls=":")
        ax1.text(min(gave_up), 0.02, " exact CASCI\n infeasible here ",
                 transform=ax1.get_xaxis_transform(), fontsize=7.5, color=INK2,
                 va="bottom", ha="right")
    ax1.set_yscale("log")
    ax1.set_xlabel("active orbitals")
    ax1.set_ylabel("wall time (s)")
    ax1.set_title("Wall clock vs active orbitals")
    ax1.legend(fontsize=7.5, loc="upper left")
    _tidy(ax1)

    cas_dims = [r["cas_dim"] for r in rows]
    ax2.plot(cas_dims, solve, "s-", color=C1, lw=2.2, ms=7, label="SQD solve")
    ax2.plot([r["cas_dim"] for r in rows if "wall_s" in r["casci"]], cy, "o-",
             color=C2, lw=2.2, ms=8, label="exact CASCI")
    # subspace size is a different measure, so it is annotated rather than plotted on
    # the same axis - one axis, one unit.
    # subspace labels go above the SQD markers with an opaque backing: below them they
    # landed on the CASCI curve where the two cross.
    for r, y in zip(rows, solve):
        ax2.annotate(f"{r['sqd_subspace_dim']:,}", (r["cas_dim"], y),
                     textcoords="offset points", xytext=(0, 11), ha="center",
                     fontsize=6.5, color=INK2,
                     bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none",
                               alpha=0.85))
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel("full CAS dimension (determinants)")
    ax2.set_ylabel("wall time (s)")
    ax2.set_title("Wall clock vs full CAS dimension")
    ax2.legend(fontsize=7.5, loc="upper left")
    _tidy(ax2)

    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def fig_orbital_opt(s2b, path):
    """Energy at *fixed* subspace dimension, before vs after orbital optimization."""
    e0 = s2b["e_casci_singlet"]
    rows = s2b["rows"]
    dim = [r["subspace_dim"] for r in rows]
    before = [(r["energy_before_oo"] - e0) * 1e3 for r in rows]
    after = [(r["energy_after_oo"] - e0) * 1e3 for r in rows]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.4, 3.8))

    ax1.plot(dim, before, "o-", color=C1, lw=2, ms=7, label="canonical AVAS basis")
    ax1.plot(dim, after, "s-", color=C2, lw=2, ms=7, label="after orbital optimization")
    ax1.axhspan(0, CHEMICAL_ACCURACY_HA * 1e3, color=C3, alpha=0.10, zorder=0)
    ax1.text(0.02, CHEMICAL_ACCURACY_HA * 1e3, " chemical accuracy",
             transform=ax1.get_yaxis_transform(), va="bottom", fontsize=7.5, color=INK2)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel("subspace dimension (frozen)")
    ax1.set_ylabel("energy error vs CASCI (mHa)")
    ax1.set_title("Energy error at fixed subspace", fontsize=10.5)
    ax1.legend(fontsize=8)
    _tidy(ax1)

    colors = [C1, C2, C3, C4, C7]
    for color, r in zip(colors, rows):
        trace = [(e - e0) * 1e3 for e in r["oo_trace"]]
        ax2.plot(range(len(trace)), trace, "o-", color=color, lw=2, ms=5,
                 label=f"dim {r['subspace_dim']}")
    ax2.set_yscale("log")
    ax2.set_xlabel("orbital-optimization iteration")
    ax2.set_ylabel("energy error (mHa)")
    ax2.set_title("Energy error vs optimisation iteration", fontsize=10.5)
    ax2.legend(fontsize=8)
    _tidy(ax2)

    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def write_markdown(s1, s2, path, figure_names, s2b=None):
    e0 = s2["e_casci_singlet"]
    fit = s1["heisenberg"]
    best = best_run(s2)
    floor = s2.get("solver_floor", {})
    L = []
    a = L.append

    a("# SQD for a [2Fe-2S(SMe)4]2- model cluster - results\n")
    a("Generated by `src/stage3_report.py`. All energies in hartree; all errors")
    a("relative to the CASCI (= full CI in the active space) singlet.\n")

    a("## System\n")
    a("| quantity | value |")
    a("|---|---|")
    a(f"| basis | {s1['basis']} |")
    a(f"| atoms / AOs / electrons | 24 / {s1['nao']} / {sum(s1['nelec_total'])} |")
    a(f"| active space | ({sum(s1['nelec_active'])}e, {s1['norb']}o) Fe 3d via AVAS "
      f"(threshold {s1['avas_threshold']}) |")
    a(f"| qubits | {s1['n_qubits']} |")
    a(f"| CAS dimension | {s1['cas_dim']} determinants |")
    a(f"| RHF | {s1['e_rhf']:.8f} |")
    a(f"| CCSD (frozen core) | {s1['e_ccsd']:.8f} "
      f"{'' if s1['ccsd_converged'] else '**(not converged)**'} |")
    a(f"| **CASCI singlet (ground truth)** | **{e0:.8f}** |")
    a("")

    a("## Baselines\n")
    a("Three tiers, in the order they matter: what we ran classically, what the")
    a("literature reports (classical *and* quantum), and what we ran on the quantum")
    a("side. Energies in hartree; `J` in cm^-1 with the Fe-S convention")
    a("`H = 2 J S1.S2`, `E(S) = J S(S+1)`, `J > 0` antiferromagnetic.\n")

    a("### Tier 1 - classical baselines, this work\n")
    a("| method | energy (Ha) | error vs CASCI | notes |")
    a("|---|---|---|---|")
    a(f"| RHF | {s1['e_rhf']:.6f} | {(s1['e_rhf'] - e0) * 1e3:+.1f} mHa | "
      f"closed-shell singlet, needed a second-order solver to converge |")
    if s1.get("e_ccsd") is not None:
        a(f"| CCSD (frozen core) | {s1['e_ccsd']:.6f} | "
          f"{(s1['e_ccsd'] - e0) * 1e3:+.1f} mHa | "
          f"**did not converge** - the diagnostic that this system is outside "
          f"single-reference range |")
    a(f"| **CASCI ({sum(s1['nelec_active'])}e,{s1['norb']}o)** | **{e0:.6f}** | "
      f"**0 (exact in-space reference)** | full CI in the active space; "
      f"{s2['cas_dim']} determinants |")
    o_mid = s2["oracle"][len(s2["oracle"]) // 2]
    a(f"| selected-CI, largest-|c| dets | {o_mid['energy']:.6f} | "
      f"{(o_mid['energy'] - e0) * 1e3:+.1f} mHa | at subspace dim "
      f"{o_mid['subspace_dim']}; the classical analogue of SQD |")
    if fit.get("ok"):
        a(f"| CASCI Heisenberg fit | - | - | **J = {fit['J_cm1']:.1f} cm^-1** |")
    a("")

    a("### Tier 2 - published reference values for this cluster\n")
    a("| method | J / accuracy | source | notes |")
    a("|---|---|---|---|")
    for name, src, val, note in LITERATURE["classical"]:
        a(f"| {name} | {val} | {src} | {note} |")
    for name, src, val, note in LITERATURE["quantum"]:
        a(f"| *(quantum)* {name} | {val} | {src} | {note} |")
    a("")

    a("### Tier 3 - quantum baselines, this work\n")
    bi = s2.get("backend_info", {})
    a(f"Ansatz: LUCJ from frozen-core CCSD `t1`/`t2`, `n_reps="
      f"{s2['config'].get('n_reps')}`, "
      f"`optimize={s2['config'].get('optimize_ansatz')}`. "
      f"Sampled on `{bi.get('backend', 'ffsim')}` "
      f"({bi.get('noise', 'noiseless')}), {s2['config'].get('shots')} shots. "
      f"Circuit depth {s2['circuit_depth']}, "
      f"{s2['circuit_ops'].get('cx', 0)} two-qubit gates.\n")
    a("| run | energy (Ha) | error vs CASCI | subspace dim | <S^2> |")
    a("|---|---|---|---|---|")
    a(f"| LUCJ ansatz expectation value (no SQD) | {s2['ansatz']['e_ansatz']:.6f} | "
      f"{(s2['ansatz']['e_ansatz'] - e0) * 1e3:+.1f} mHa | - | - |")
    a(f"| **SQD, best** | **{best['energy']:.6f}** | "
      f"**{(best['energy'] - e0) * 1e3:+.4f} mHa** | {best['subspace_dim']} | "
      f"{best['spin_sq']:.4f} |")
    u = s2["uniform_control"]
    a(f"| SQD on uniform-random bitstrings (control) | {u['energy']:.6f} | "
      f"{(u['energy'] - e0) * 1e3:+.4f} mHa | {u['subspace_dim']} | "
      f"{u['spin_sq']:.4f} |")
    if s2b:
        bb = min(s2b["rows"], key=lambda r: r["energy_after_oo"])
        a(f"| SQD + orbital optimization (frozen subspace) | "
          f"{bb['energy_after_oo']:.6f} | "
          f"{(bb['energy_after_oo'] - e0) * 1e3:+.1f} mHa | "
          f"{bb['subspace_dim']} | - |")
    a("")
    a("**Verdict.** On energy we are *exact* against the in-space reference, which is a")
    a("stronger result than the published SQD run on this cluster (tens of mHa vs HCI) -")
    a("but only because our active space is small enough to diagonalize exactly, so this")
    a("is a validation of the pipeline, not a demonstration of quantum advantage.")
    if fit.get("ok"):
        ref = fit["reference_J_cm1"]["DMRG (30e,20o), Sharma/Chan 2014"]
        a(f"On the metric the field actually reports, **J = {fit['J_cm1']:.1f} cm^-1 "
          f"against DMRG's {ref:.0f} and experiment's 148 +/- 16** - low by roughly a")
        a("factor of eight. The cause is the active space and basis, not the solver:")
        a("Fe 3d only means no bridging-sulfur 3p superexchange pathway, and STO-3G is a")
        a("minimal basis where the literature uses TZP-DKH.")
    a("")

    a("## Headline\n")
    a(f"- Best SQD energy: **{best['energy']:.8f}** "
      f"(**{(best['energy'] - e0) * 1e3:+.3f} mHa**) at subspace dimension "
      f"{best['subspace_dim']} "
      f"({best['subspace_dim'] / s2['cas_dim'] * 100:.1f}% of the CAS space), "
      f"⟨S²⟩ = {best['spin_sq']:.4f}.")
    if floor:
        a(f"- Diagonalizing the *entire* CAS space through the same selected-CI solver "
          f"leaves {floor['error_mha']:+.3f} mHa, so that is the floor for every number "
          f"here - not a property of the sampling.")
    a(f"- The LUCJ ansatz itself sits {(s2['ansatz']['e_ansatz'] - e0) * 1e3:+.1f} mHa "
      f"above CASCI with a participation ratio of "
      f"{s2['ansatz']['participation_ratio']:.2f}; SQD's accuracy comes from the "
      f"classical diagonalization, not from the ansatz energy.")
    a(f"- A uniform-random-bitstring control reaches "
      f"{(s2['uniform_control']['energy'] - e0) * 1e3:+.3f} mHa at dimension "
      f"{s2['uniform_control']['subspace_dim']}. At this active-space size the "
      f"circuit's distribution carries no advantage over random configurations.")
    a("")

    a("## Energy vs subspace dimension\n")
    a(f"![convergence]({figure_names['convergence']})\n")
    a("| samples/batch | subspace dim | % of CAS | E (Ha) | error (mHa) | ⟨S²⟩ |")
    a("|---|---|---|---|---|---|")
    for r in s2["dim_sweep"]:
        a(f"| {r['samples_per_batch']} | {r['subspace_dim']} | "
          f"{r['subspace_dim'] / s2['cas_dim'] * 100:.1f}% | {r['energy']:.8f} | "
          f"{(r['energy'] - e0) * 1e3:+.4f} | {r['spin_sq']:.4f} |")
    a("")
    a("Selected-CI reference: keep the largest-|c| CASCI determinants and diagonalize in")
    a("the product closure of their spin-up and spin-down strings. Note this is a")
    a("*reference strategy*, not an upper bound - the product closure drags in many")
    a("low-weight configurations, and at small dimensions SQD's sampled subspace actually")
    a("beats it:\n")
    a("| dets kept | subspace dim | error (mHa) |")
    a("|---|---|---|")
    for r in s2["oracle"]:
        a(f"| {r['n_dets_kept']} | {r['subspace_dim']} | {(r['energy'] - e0) * 1e3:+.4f} |")
    a("")

    a("## Sampling efficiency vs noise\n")
    a(f"![noise]({figure_names['noise']})\n")
    a("| depolarizing rate | valid shots (correct N and Sz) | unique bitstrings | "
      "subspace dim | error (mHa) | ⟨S²⟩ |")
    a("|---|---|---|---|---|---|")
    for r in s2["noise_sweep"]:
        e = r["efficiency"]
        a(f"| {r['depol']:g} | {e['frac_correct_particle_number_and_sz'] * 100:.2f}% | "
          f"{e['unique_bitstrings']} | {r['subspace_dim']} | "
          f"{(r['energy'] - e0) * 1e3:+.4f} | {r['spin_sq']:.4f} |")
    a("")

    a("## Configuration-recovery convergence\n")
    a(f"![iterations]({figure_names['iterations']})\n")

    a("## Orbital occupancies\n")
    a(f"![occupancies]({figure_names['occupancies']})\n")
    exact = np.array(s2["exact_occupancies"]["occ_a"]) + np.array(
        s2["exact_occupancies"]["occ_b"])
    sqd = np.array(best["occ_a"]) + np.array(best["occ_b"])
    a(f"Maximum per-orbital deviation from CASCI: **{np.abs(sqd - exact).max():.4f}** "
      f"electrons; the exact occupancies are all close to 1.0, i.e. a half-filled, "
      f"open-shell d manifold on both irons.\n")

    if s2b:
        a("## Orbital optimization at fixed subspace\n")
        a("The subspace is frozen to the configurations SQD found and only the orbital")
        a("basis is varied, so any energy drop here is genuine compression rather than a")
        a("larger subspace (port of the package's `optimize_orbitals` guide).\n")
        a(f"![orbital optimization]({figure_names['orbital_opt']})\n")
        a("| subspace dim | % of CAS | error before OO (mHa) | error after OO (mHa) | gain |")
        a("|---|---|---|---|---|")
        for r in s2b["rows"]:
            eb = (r["energy_before_oo"] - e0) * 1e3
            ea = (r["energy_after_oo"] - e0) * 1e3
            gain = f"{eb / ea:.1f}x" if abs(ea) > 1e-9 else "-"
            a(f"| {r['subspace_dim']} | {r['subspace_dim'] / s2b['cas_dim'] * 100:.1f}% "
              f"| {eb:+.3f} | {ea:+.3f} | {gain} |")
        a("")
        bestb = min(s2b["rows"], key=lambda r: r["energy_after_oo"])
        eb = (bestb["energy_before_oo"] - e0) * 1e3
        ea = (bestb["energy_after_oo"] - e0) * 1e3
        a(f"Best after optimization: **{ea:+.3f} mHa** at subspace dimension "
          f"{bestb['subspace_dim']} "
          f"({bestb['subspace_dim'] / s2b['cas_dim'] * 100:.1f}% of the CAS space), "
          f"from {eb:+.3f} mHa at the same dimension in the canonical AVAS basis.\n")
        pr = s2["exact_wavefunction"]["participation_ratio"]
        n99 = s2["exact_wavefunction"]["dets_for_weight"]["0.99"]
        a(f"**Orbital optimization does not rescue the compression here.** It recovers "
          f"{eb - ea:.1f} mHa of a {eb:.0f} mHa gap ({(eb - ea) / eb * 100:.0f}%), and "
          f"the traces flatten after a few iterations. That is consistent with the")
        a(f"structure of the exact state rather than a bad choice of basis: in this basis")
        a(f"the CASCI singlet has a participation ratio of {pr:.0f} and needs {n99} "
          f"determinants for 99% of its weight, out of {s2['cas_dim']}. A wavefunction")
        a("that flat has no compact product subspace to find, so neither a better basis")
        a("nor more shots helps much - only covering most of the space does.\n")

    a("## Physics: the spin ladder and the exchange coupling\n")
    a(f"![spin ladder]({figure_names['spin_ladder']})\n")
    a("| S | ⟨S²⟩ | E (Ha) | above singlet (mHa) |")
    a("|---|---|---|---|")
    emin = min(st["e_total"] for st in s1["casci_states"])
    for st in s1["casci_states"]:
        a(f"| {st['S']:.0f} | {st['spin_sq']:.4f} | {st['e_total']:.8f} | "
          f"{(st['e_total'] - emin) * 1e3:.4f} |")
    a("")
    if fit.get("ok"):
        a(f"The six roots follow the Landé interval rule to "
          f"{fit['max_residual_hartree']:.1e} Ha, giving **J = {fit['J_cm1']:.1f} cm⁻¹** "
          f"(convention {fit['convention']}). Oxidized synthetic and protein [2Fe-2S] "
          f"centres sit near -150 to -400 cm⁻¹, so an Fe-3d-only active space in this "
          f"basis underestimates the antiferromagnetic coupling by roughly an order of "
          f"magnitude. The bridging-sulfur 3p superexchange pathways are missing from "
          f"the active space, and no amount of extra sampling can recover them.\n")

    a("## What this says about SQD here\n")
    a("1. **The pipeline reproduces exact diagonalization.** Post-selection plus")
    a("   self-consistent configuration recovery drives the energy to the solver floor")
    a("   and the occupancies to the CASCI values, with ⟨S²⟩ driven to ~0.")
    a("2. **It does so by covering essentially the whole CAS space, not by compressing")
    a("   it.** The subspace SQD builds is a Cartesian product of the sampled spin-up")
    a("   and spin-down strings, and in the canonical AVAS orbital basis the exact")
    a(f"   singlet has a participation ratio of "
      f"{s2['exact_wavefunction']['participation_ratio']:.0f} with "
      f"{s2['exact_wavefunction']['dets_for_weight']['0.99']} determinants needed for 99%")
    a("   of the weight. Any sampler with modest diversity therefore reconstructs the")
    a("   full space, which is why the uniform-random control matches the circuit.")
    a("3. **The 1.6 mHa target is the wrong yardstick for this system.** The entire")
    a(f"   spin ladder spans {(max(st['e_total'] for st in s1['casci_states']) - emin) * 1e3:.2f} mHa,")
    a("   so chemical accuracy does not distinguish a singlet from a septet here. The")
    a("   meaningful target is the ladder spacing itself, which is what J depends on.")
    a("4. **The binding constraint is the active space, not the solver.** Getting J")
    a("   right needs the bridging-S 3p orbitals; that is the direction a larger run")
    a("   should go, and it is also where the CAS space becomes too large to")
    a("   diagonalize exactly - i.e. where SQD would start to earn its keep.")
    a("")
    a("### Caveats\n")
    a(f"- The frozen-core CCSD that seeds the LUCJ amplitudes **did not converge** "
      f"(|t2|max is large) - expected for a half-filled 3d shell, and the reason the")
    a("  ansatz energy is poor. The amplitudes are still usable as a sampling")
    a("  distribution, which is all SQD asks of them.")
    a("- Hardware noise is emulated with ffsim's global depolarizing channel rather")
    a("  than a device noise model, so the valid-shot fractions are indicative only.")
    a("- The geometry is an idealized model cluster, not a relaxed or crystal structure.")
    a("")

    path.write_text("\n".join(L) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="sto-3g_avas0.5")
    ap.add_argument("--suffix", default="")
    args = ap.parse_args()

    s1 = load_json(RESULTS / f"stage1_{args.tag}.json")
    s2 = load_json(RESULTS / f"stage2_{args.tag}{args.suffix}.json")
    FIGS.mkdir(parents=True, exist_ok=True)

    names = {
        "convergence": f"figures/convergence_{args.tag}{args.suffix}.png",
        "noise": f"figures/noise_{args.tag}{args.suffix}.png",
        "iterations": f"figures/iterations_{args.tag}{args.suffix}.png",
        "occupancies": f"figures/occupancies_{args.tag}{args.suffix}.png",
        "spin_ladder": f"figures/spin_ladder_{args.tag}.png",
        "orbital_opt": f"figures/orbital_opt_{args.tag}.png",
    }
    fig_convergence(s2, RESULTS / names["convergence"])
    fig_noise(s2, RESULTS / names["noise"])
    fig_iterations(s2, RESULTS / names["iterations"])
    fig_occupancies(s2, RESULTS / names["occupancies"])
    fig_spin_ladder(s1, RESULTS / names["spin_ladder"])

    s5_path = RESULTS / "stage5_scaling.json"
    if s5_path.exists():
        names["scaling"] = "figures/scaling_casci_vs_sqd.png"
        fig_scaling(load_json(s5_path), RESULTS / names["scaling"])

    s2b_path = RESULTS / f"stage2b_{args.tag}.json"
    s2b = load_json(s2b_path) if s2b_path.exists() else None
    if s2b:
        fig_orbital_opt(s2b, RESULTS / names["orbital_opt"])
    else:
        names.pop("orbital_opt")

    out = RESULTS / f"RESULTS_{args.tag}{args.suffix}.md"
    write_markdown(s1, s2, out, names, s2b=s2b)
    print("wrote", out)
    for n in names.values():
        print("wrote", RESULTS / n)


if __name__ == "__main__":
    main()

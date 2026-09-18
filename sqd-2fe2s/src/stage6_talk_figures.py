"""Stage 6 - the three figures for a talk.

Deliberately separate from stage3 (which documents everything). These are built to be
understood in ten seconds by someone who has never heard of an active space.
"""
from __future__ import annotations
import json, pathlib, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from common import RESULTS  # noqa: E402

C1, C2, C3, C4, C7 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#4a3aa7"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8984"
plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 150, "font.size": 10,
    "axes.edgecolor": INK3, "axes.labelcolor": INK, "axes.titlesize": 11.5,
    "axes.titleweight": "semibold", "axes.grid": True, "grid.color": "#e4e3df",
    "grid.linewidth": 0.6, "xtick.color": INK2, "ytick.color": INK2,
    "legend.frameon": False, "axes.spines.top": False, "axes.spines.right": False})
FIGS = RESULTS / "figures"; FIGS.mkdir(parents=True, exist_ok=True)


def fig_runtime():
    """Runtime vs problem size: exact diagonalisation against SQD's solve.

    Uses the measured nested-slice study.  Both curves are wall clock on the same
    machine, one core, so the comparison at a given size is fair.
    """
    s5 = json.load(open(RESULTS / "stage5_scaling.json"))
    rows = s5["rows"]
    dims = [r["cas_dim"] for r in rows]
    solve = [r["sqd_solve_s"] for r in rows]
    cx = [r["cas_dim"] for r in rows if "wall_s" in r["casci"]]
    cy = [r["casci"]["wall_s"] for r in rows if "wall_s" in r["casci"]]

    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    ax.plot(cx, cy, "o-", color=C2, lw=2.6, ms=9,
            label="Exact diagonalisation (CASCI)", zorder=5)
    ax.plot(dims, solve, "s-", color=C1, lw=2.6, ms=8,
            label="SQD — classical solve", zorder=6)

    wall = [r["cas_dim"] for r in rows if "wall_s" not in r["casci"]]
    if wall:
        ax.axvline(min(wall), color=INK3, lw=1.4, ls=":")
        ax.annotate("exact runs out of\nmemory here",
                    xy=(min(wall), max(cy)), xytext=(-12, -10),
                    textcoords="offset points", ha="right", va="top",
                    fontsize=8.5, color=INK2)

    for sec, lab in [(60, "1 min"), (3600, "1 hour")]:
        ax.axhline(sec, color=INK3, lw=0.8, ls="--", alpha=0.6)
        ax.annotate(lab, xy=(dims[0], sec), xytext=(2, 3), ha="left",
                    textcoords="offset points", fontsize=8, color=INK2)

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("problem size (number of electron arrangements)")
    ax.set_ylabel("time (seconds)")
    ax.set_title("Wall clock vs problem size: exact diagonalisation and the SQD solve\n"
                 "nested slices of the [2Fe-2S] Hamiltonian, one core")
    ax.legend(loc="lower right", fontsize=9)

    # Second axis in qubits, because that is the unit people ask the question in.
    qb = [r["n_qubits"] for r in rows]
    top = ax.secondary_xaxis("top")
    top.set_xticks(dims)
    top.set_xticklabels([f"{q}q" for q in qb], fontsize=8.5, color=INK2)
    top.tick_params(length=2)
    ax.set_axisbelow(True)
    fig.tight_layout(); fig.savefig(FIGS / "talk_runtime.png", bbox_inches="tight")
    plt.close(fig); print("wrote figures/talk_runtime.png")


def fig_energy_vs_iteration():
    """Energy against iteration: the self-consistent loop closing on the exact answer."""
    d = json.load(open(RESULTS / "timing_22e16o_clean.json"))
    h = d["history"]
    it = [r["iter"] for r in h]
    err = [r["err_mha"] for r in h]
    dim = [r["dim"] for r in h]
    E_EXACT = -5013.72724223
    energy = [E_EXACT + e / 1e3 for e in err]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.0, 6.2), sharex=True,
                                   gridspec_kw={"height_ratios": [2, 1]})

    ax1.axhline(E_EXACT, color=C2, lw=2.2, ls="--",
                label=f"exact answer  {E_EXACT:.5f} Ha")
    ax1.plot(it, energy, "o-", color=C1, lw=2.4, ms=8, label="SQD", zorder=5)
    ax1.set_ylabel("energy (Hartree)")
    ax1.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.2f}"))
    ax1.set_title("Energy vs configuration-recovery iteration\n"
                  "(22e,16o), 32 qubits, 19,079,424 determinants")
    ax1.legend(loc="upper right", fontsize=9)
    for i in (0, len(it) // 2, len(it) - 1):
        ax1.annotate(f"{err[i]:+.0f} mHa", (it[i], energy[i]),
                     textcoords="offset points", xytext=(6, 8), fontsize=8, color=INK2)
    ax1.set_axisbelow(True)

    ax2.plot(it, dim, "d-", color=C7, lw=2.0, ms=7)
    ax2.set_xlabel("iteration of the self-consistent loop")
    ax2.set_ylabel("configurations kept")
    ax2.set_title("Subspace dimension per iteration", fontsize=10)
    ax2.set_xticks(it); ax2.set_axisbelow(True)

    fig.tight_layout(); fig.savefig(FIGS / "talk_energy_vs_iteration.png",
                                    bbox_inches="tight")
    plt.close(fig); print("wrote figures/talk_energy_vs_iteration.png")


def fig_j_vs_active_space():
    """The physics payoff: J converging toward experiment as the model improves."""
    spaces = ["Fe 3d only\n(20 qubits)", "+ bridging S\n(32 qubits)",
              "+ Fe 4s\n(36 qubits)"]
    J = [28.8, 46.7, 82.0]
    x = np.arange(len(J))

    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    ax.axhspan(132, 164, color=C3, alpha=0.16, zorder=0)
    ax.axhline(148, color=C3, lw=2.4)
    ax.axhline(236, color=INK3, lw=1.6, ls="--")
    ax.bar(x, J, width=0.55, color=C1, edgecolor="white", linewidth=1.2, zorder=4)
    for xi, j in zip(x, J):
        ax.annotate(f"{j:.0f}", (xi, j), textcoords="offset points", xytext=(0, 5),
                    ha="center", fontsize=10.5, color=INK, fontweight="semibold")

    # The reference lines are labelled inline, in the margin to the right of the bars:
    # a legend box here sat on top of the 236 line whatever corner it was put in.
    ax.set_xlim(-0.62, len(J) - 0.08)
    for y, txt, col in [(236, "best published\ncalculation  236", INK2),
                        (148, "experiment\n148 ± 16", C3)]:
        ax.annotate(txt, xy=(len(J) - 0.12, y), xytext=(-2, 0),
                    textcoords="offset points", ha="right", va="center",
                    fontsize=9, color=col, fontweight="semibold")
    ax.annotate("this work", xy=(0, J[0]), xytext=(-46, 10),
                textcoords="offset points", ha="center", fontsize=9.5, color=C1,
                fontweight="semibold")

    ax.set_xticks(x); ax.set_xticklabels(spaces, fontsize=9)
    ax.set_ylabel("magnetic coupling J (cm⁻¹)")
    ax.set_title("Exchange coupling J by active space")
    ax.set_ylim(0, 270); ax.set_axisbelow(True)
    fig.tight_layout(); fig.savefig(FIGS / "talk_J_vs_active_space.png",
                                    bbox_inches="tight")
    plt.close(fig); print("wrote figures/talk_J_vs_active_space.png")


if __name__ == "__main__":
    fig_runtime(); fig_energy_vs_iteration(); fig_j_vs_active_space()

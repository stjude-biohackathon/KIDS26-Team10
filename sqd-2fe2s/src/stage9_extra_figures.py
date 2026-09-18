"""Stage 9 - three figures for results we had measured but never plotted.

All from existing JSON, no recomputation:

    python3 src/stage9_extra_figures.py

  coverage_22e16o.png   energy error vs fraction of the space the subspace covers
  cost_casci_vs_sqd.png wall clock and per-iteration error, exact vs SQD at 32 qubits
  shots_null.png        unique bitstrings and energy error vs shot count
"""
from __future__ import annotations

import json
import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from common import RESULTS  # noqa: E402

C1, C2, C3, C4, C7 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#4a3aa7"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8984"
plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 150, "font.size": 10,
    "axes.edgecolor": INK3, "axes.labelcolor": INK, "axes.titlesize": 11.5,
    "axes.titleweight": "semibold", "axes.grid": True, "grid.color": "#e4e3df",
    "grid.linewidth": 0.6, "xtick.color": INK2, "ytick.color": INK2,
    "legend.frameon": False, "axes.spines.top": False, "axes.spines.right": False})
FIGS = RESULTS / "figures"
FIGS.mkdir(parents=True, exist_ok=True)


def fig_coverage():
    """Error against the fraction of the 19.1M-determinant space the subspace covers.

    Two series, because they are two different configurations and merging them into one
    curve would be dishonest:
      * the subspace sweep - K = 1 batch, 2 recovery iterations, 200k shots
      * the long run       - K = 3 batches, 10 recovery iterations, same shots
    The point they make together is the same: on this molecule the error only becomes
    interesting once the subspace is a large fraction of the whole space.
    """
    curve = json.load(open(RESULTS / "stage2_sto-3g_fe3d+brs3p_curve.json"))
    long = json.load(open(RESULTS / "timing_22e16o_clean.json"))
    pct = [r["pct_cas"] for r in curve["rows"]]
    err = [r["err_mha"] for r in curve["rows"]]

    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    ax.plot(pct, err, "o-", color=C1, lw=2.6, ms=9, zorder=5,
            label="subspace sweep (K = 1, 2 recovery iterations)")
    ax.plot([long["pct_cas"]], [long["err_mha"]], "D", color=C2, ms=13, zorder=6,
            label="best run (K = 3, 10 recovery iterations)")
    # CASCI is exactly 0 mHa, which a log axis cannot show - pin the marker to the
    # bottom of the axis and label it, rather than implying a small finite error.
    ax.plot([100], [1.15], "*", color=C3, ms=24, zorder=6, clip_on=False,
            label="100% of the space = exact diagonalisation, 0 mHa")

    ax.axhline(1.6, color=C7, lw=1.8, ls="--", zorder=4)
    ax.text(0.115, 1.12, "chemical accuracy 1.6 mHa", color=C7, fontsize=9)
    ax.annotate(f"{long['pct_cas']:.0f}% of the space\n+{long['err_mha']:.1f} mHa",
                (long["pct_cas"], long["err_mha"]), textcoords="offset points",
                xytext=(-4, -54), fontsize=10.5, color=C2, fontweight="semibold",
                ha="center", arrowprops=dict(arrowstyle="->", color=C2, lw=1.3))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("subspace as a fraction of the full space (%, log scale)")
    ax.set_ylabel("error vs exact diagonalisation (mHa, log scale)")
    ax.set_title("Energy error vs subspace coverage\n"
                 "(22e,16o), 32 qubits, 19,079,424 determinants")
    ax.set_xlim(0.1, 220)
    ax.set_ylim(1, 3000)
    ax.legend(loc="center left", fontsize=9)
    ax.text(.985, .97,
            "participation ratio ~1,956;\n"
            "7,171 determinants hold 99%\nof the wavefunction weight",
            transform=ax.transAxes, fontsize=8.8, color=INK2, ha="right", va="top")
    fig.tight_layout()
    p = FIGS / "coverage_22e16o.png"
    fig.savefig(p)
    plt.close(fig)
    print(f"  wrote {p.name}")


def fig_cost():
    """Wall clock at 32 qubits: one exact Davidson against ten SQD recovery iterations.

    Both measured on the same otherwise-idle machine, one core, so this is a fair
    comparison at this problem size - and exact diagonalisation wins it.
    """
    t = json.load(open(RESULTS / "timing_22e16o_clean.json"))
    casci_h = t["casci_singlet_s"] / 3600
    sqd_h = t["t_sqd_s"] / 3600

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(8.6, 4.7),
                                  gridspec_kw={"width_ratios": [1, 1.35]})

    bars = ax.bar([0, 1], [casci_h, sqd_h], color=[C3, C2], width=.6, zorder=3)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["exact\n(CASCI)", "SQD\n(10 iterations)"])
    ax.set_ylabel("wall clock, one core (hours)")
    for x, v, lab in zip([0, 1], [casci_h, sqd_h],
                         [f"{casci_h:.2f} h\nexact", f"{sqd_h:.2f} h\n+{t['err_mha']:.1f} mHa"]):
        ax.annotate(lab, (x, v), textcoords="offset points", xytext=(0, 6),
                    ha="center", fontsize=10, color=INK, fontweight="semibold")
    ax.set_ylim(0, sqd_h * 1.35)
    ax.set_title("Wall clock")

    it = [h["iter"] for h in t["history"]]
    er = [h["err_mha"] for h in t["history"]]
    ax2.plot(it, er, "o-", color=C2, lw=2.6, ms=8, zorder=5)
    ax2.axhline(1.6, color=C7, lw=1.8, ls="--")
    ax2.text(10, 2.1, "1.6 mHa", color=C7, fontsize=9, ha="right")
    ax2.set_yscale("log")
    ax2.set_xlabel("configuration-recovery iteration")
    ax2.set_ylabel("error vs exact (mHa, log scale)")
    ax2.set_title("Energy error vs iteration")
    ax2.set_xticks(it)
    ax2.set_ylim(1, 900)
    ax2.annotate(f"+{er[-1]:.1f} mHa\nat the cap", (it[-1], er[-1]),
                 textcoords="offset points", xytext=(-16, 26), fontsize=9.5,
                 color=C2, fontweight="semibold", ha="center",
                 arrowprops=dict(arrowstyle="->", color=C2, lw=1.2))
    fig.tight_layout(rect=(0, 0.075, 1, 1))
    fig.text(.5, .02, "Exact diagonalisation and SQD at 32 qubits - (22e,16o), "
             "19,079,424 determinants, one core, idle machine.", ha="center",
             fontsize=9.5, color=INK, fontweight="semibold")
    p = FIGS / "cost_casci_vs_sqd.png"
    fig.savefig(p)
    plt.close(fig)
    print(f"  wrote {p.name}")


def fig_shots_null():
    """Shot count against energy: 12x more shots, and the answer got slightly worse.

    Same subspace configuration in both arms - only the shot count changed.  The reason
    is `samples_per_batch`, which caps how many configurations enter the subspace
    regardless of how good the sample pool is.
    """
    ab = json.load(open(RESULTS / "shots_ab.json"))
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(8.6, 4.6))
    lab = [f"{a['shots']:,}\nshots" for a in ab]
    x = np.arange(2)

    ax.bar(x, [a["unique"] for a in ab], color=C1, width=.55, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(lab)
    ax.set_ylabel("unique bitstrings sampled")
    ax.set_title("Unique bitstrings sampled")
    for xi, a in zip(x, ab):
        ax.annotate(f"{a['unique']:,}", (xi, a["unique"]), ha="center",
                    textcoords="offset points", xytext=(0, 5), fontsize=10,
                    color=INK, fontweight="semibold")
    ax.set_ylim(0, max(a["unique"] for a in ab) * 1.18)

    er = [a["err_mha"] for a in ab]
    ax2.bar(x, er, color=[INK3, C2], width=.55, zorder=3)
    ax2.set_xticks(x)
    ax2.set_xticklabels(lab)
    ax2.set_ylabel("error vs exact diagonalisation (mHa)")
    ax2.set_title("Energy error vs exact diagonalisation")
    for xi, e in zip(x, er):
        ax2.annotate(f"+{e:.1f} mHa", (xi, e), ha="center",
                     textcoords="offset points", xytext=(0, 5), fontsize=10.5,
                     color=INK, fontweight="semibold")
    ax2.set_ylim(0, max(er) * 1.3)
    fig.tight_layout(rect=(0, 0.10, 1, 1))
    fig.text(.5, .055, "Shot count vs sampling and accuracy - (22e,16o), 32 qubits, "
             "identical subspace configuration in both arms.", ha="center",
             fontsize=9.5, color=INK, fontweight="semibold")
    fig.text(.5, .015, "samples_per_batch caps how many configurations enter the "
             "subspace, independent of the sample pool.",
             ha="center", fontsize=9, color=INK2)
    p = FIGS / "shots_null.png"
    fig.savefig(p)
    plt.close(fig)
    print(f"  wrote {p.name}")


if __name__ == "__main__":
    fig_coverage()
    fig_cost()
    fig_shots_null()
    print("done")

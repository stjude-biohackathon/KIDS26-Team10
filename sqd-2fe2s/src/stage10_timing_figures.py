"""Stage 10 - the two figures a collaborator keeps asking for, from measured wall clock.

    python3 src/stage10_timing_figures.py

  runtime_vs_active_space.png   CASCI vs SQD wall clock, 20 -> 32 -> 36 qubits, measured
  gpu_casci_sectors.png         GPU exact-CASCI time vs sector dimension at (26e,18o)
  energy_vs_iteration_20q.png   energy in HARTREE vs recovery iteration, 20 qubits

Every number is read from a file in results/ or is listed in TRACE below with the file it
came from, so any of them can be checked. Nothing here recomputes anything.
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

# Measured wall clock on the three REAL active spaces.  Two numbers are not in a JSON
# file, so they are pinned here with their source:
#   20q CASCI 11.9 s  - full 6-root spin ladder, from the stage-1 run log quoted in
#                       docs/RUN_SUMMARY.md section 2.  Singlet alone is ~0.4 s.
#   20q SQD   84.9 s  - results/stage2_sto-3g_fe3d_opt.json, dim_sweep spb=400, the
#                       cheapest configuration that actually reaches the exact answer.
TRACE = {
    20: dict(space="(10e,10o)", dim=63504, casci_s=11.9, sqd_s=84.9,
             casci_err=0.0, sqd_err=0.0, where="laptop, 1 core"),
    32: dict(space="(22e,16o)", dim=19079424, casci_s=None, sqd_s=None,
             casci_err=0.0, sqd_err=None, where="idle laptop, 1 core"),
    36: dict(space="(26e,18o)", dim=73410624, casci_s=None, sqd_s=None,
             casci_err=0.0, sqd_err=None, where="A100 GPU"),
}


def load():
    t32 = json.load(open(RESULTS / "timing_22e16o_clean.json"))
    TRACE[32]["casci_s"] = t32["casci_singlet_s"]
    TRACE[32]["sqd_s"] = t32["t_sqd_s"]
    TRACE[32]["sqd_err"] = t32["err_mha"]
    g36 = json.load(open(RESULTS / "stage1_26e18o_gpu_ladder.json"))
    # four exact spin sectors; there is no single-sector singlet time because Sz=0 and
    # Sz=1 never fitted in memory, so the honest number is the sum of what did run.
    TRACE[36]["casci_s"] = sum(s["s"] for s in g36["sectors"])
    TRACE[36]["sectors"] = g36["sectors"]
    return t32, g36


def fig_runtime_vs_active_space():
    """Measured wall clock against active-space size, for the three spaces we ran.

    This is the plot to use when someone asks "how does runtime grow with active space
    size" - unlike the nested-slice study in stage 5, every point here is a real
    calculation on the real Hamiltonian, so the energies attached to them mean something.

    The 36-qubit CASCI point ran on a GPU and is drawn hollow: it is NOT comparable to
    the two CPU points, and pretending otherwise would be the easiest mistake to make
    with this figure.
    """
    qs = [20, 32, 36]
    dims = [TRACE[q]["dim"] for q in qs]
    casci = [TRACE[q]["casci_s"] / 3600 for q in qs]
    sqd = [TRACE[q]["sqd_s"] / 3600 if TRACE[q]["sqd_s"] else np.nan for q in qs]

    fig, ax = plt.subplots(figsize=(7.6, 5.0))
    ax.plot(dims[:2], casci[:2], "o-", color=C3, lw=2.8, ms=11, zorder=5,
            label="exact diagonalisation (CASCI), CPU")
    ax.plot([dims[2]], [casci[2]], "o", color=C3, ms=13, mfc="white", mew=2.6,
            zorder=5, label="exact diagonalisation, GPU — not CPU-comparable")
    ax.plot(dims[:2], sqd[:2], "s-", color=C2, lw=2.8, ms=10, zorder=5,
            label="SQD, full self-consistent loop, CPU")

    def hrs(h):
        """Hours are unreadable for the small runs; switch units rather than print 0.02 h."""
        return f"{h * 3600:.0f} s" if h < 0.1 else f"{h:.2f} h"

    # 20q sits on the axis floor, and 36q has the memory note below it, so both of
    # those labels go above their marker instead of below.
    OFFSET = {20: (0, 14), 32: (0, -42), 36: (-2, 16)}
    for q, d, c in zip(qs, dims, casci):
        ax.annotate(f"{q}q  {TRACE[q]['space']}\n{hrs(c)}",
                    (d, c), textcoords="offset points", xytext=OFFSET[q],
                    ha="center", fontsize=9, color=C3, fontweight="semibold")
    for q, d, s in zip(qs[:2], dims[:2], sqd[:2]):
        err = TRACE[q]["sqd_err"]
        lab = f"{hrs(s)}\n" + ("exact" if err == 0 else f"+{err:.1f} mHa")
        ax.annotate(lab, (d, s), textcoords="offset points", xytext=(0, 12),
                    ha="center", fontsize=9.5, color=C2, fontweight="semibold")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("active-space size (determinants, log scale)")
    ax.set_ylabel("wall clock (hours, log scale)")
    ax.set_title("Wall clock vs active-space size\n"
                 "20, 32 and 36 qubits on the [2Fe-2S] Hamiltonian")
    # Memory, not cores, is what actually stopped each rung - mark where.
    ax.annotate("17.6 GB working set:\ntoo big for 16 GB of laptop RAM,\n"
                "ran on a GPU. 2 of 6 spin\nsectors still out of memory",
                xy=(dims[2], casci[2]), xytext=(-6, -96), textcoords="offset points",
                ha="center", fontsize=8.6, color=C2, fontweight="semibold",
                arrowprops=dict(arrowstyle="->", color=C2, lw=1.2))
    # fixed axes-fraction placement: offsetting from the 32q marker put this note on
    # top of the SQD line whichever direction it was nudged
    ax.text(.30, .60, "153 MB per CI vector,\n~30 kept by the Davidson",
            transform=ax.transAxes, ha="center", fontsize=8.6, color=INK2)
    ax.set_xlim(2e4, 6e8)
    ax.set_ylim(2e-3, 40)
    ax.legend(loc="upper left", fontsize=9)
    fig.tight_layout()
    p = FIGS / "runtime_vs_active_space.png"
    fig.savefig(p)
    plt.close(fig)
    print(f"  wrote {p.name}")


def fig_gpu_sectors():
    """Exact-CASCI time against sector dimension on the GPU, at (26e,18o).

    Four spin sectors, 43,758 to 25,968,384 determinants - a clean small-to-large
    runtime series inside one active space, on one device.  The two missing sectors
    (Sz = 1 and Sz = 0, 56.8M and 73.4M) ran out of memory, which is why J at this
    space is a four-point fit.
    """
    sec = TRACE[36]["sectors"]
    d = np.array([s["dim"] for s in sec], float)
    t = np.array([s["s"] for s in sec], float)
    sz = [s["Sz"] for s in sec]

    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    ax.plot(d, t, "o-", color=C7, lw=2.8, ms=11, zorder=5, label="measured, A100")
    # No extrapolation line: four points over three decades do not pin a power law well
    # enough here (the local slope steepens 0.50 -> 1.28 -> 2.01), and a fitted curve
    # would imply more than the data supports.  The two sectors that never ran are
    # marked on the x-axis instead.
    # Sz=5 is at the left edge, so its label goes right of the marker; the rest go left,
    # clear of the rising curve.
    for di, ti, s in zip(d, t, sz):
        dx = 34 if s == max(sz) else -34
        ax.annotate(f"Sz={s}\n{ti:,.0f} s", (di, ti), textcoords="offset points",
                    xytext=(dx, 2), ha="center", fontsize=9, color=C7,
                    fontweight="semibold")
    for dd in (56.8e6, 73.4e6):
        ax.axvline(dd, color=C2, lw=1.6, ls=":", zorder=3)
    ax.axvspan(50e6, 3e8, color=C2, alpha=.07, zorder=1)
    # opaque backing: this note sits over the two dotted out-of-memory markers
    ax.text(2.9e8, 4.6, "Sz = 1 (56.8M) and Sz = 0 (73.4M)\nran out of memory - which is\n"
            "why J at this space is a four-point fit",
            fontsize=8.8, color=C2, fontweight="semibold", ha="right", va="bottom",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.92))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("sector dimension (determinants, log scale)")
    ax.set_ylabel("exact-CASCI wall clock (seconds, log scale)")
    ax.set_title("Exact-CASCI wall clock vs spin-sector dimension\n"
                 "(26e,18o), 36 qubits, A100 GPU")
    ax.set_xlim(2e4, 5e8)
    ax.set_ylim(3, 3e4)
    ax.legend(loc="upper left", fontsize=9)
    fig.tight_layout()
    p = FIGS / "gpu_casci_sectors.png"
    fig.savefig(p)
    plt.close(fig)
    print(f"  wrote {p.name}")


def fig_energy_vs_iteration_20q():
    """Energy in HARTREE against configuration-recovery iteration, 20 qubits.

    Six subspace sizes on one axis, against the exact CASCI singlet.  This is the
    20-qubit companion to talk_energy_vs_iteration.png (which is the 32-qubit run):
    here the curves actually land on the exact answer, which the 32-qubit one never does.
    """
    d = json.load(open(RESULTS / "stage2_sto-3g_fe3d_opt.json"))
    exact = d["e_casci_singlet"]
    rows = sorted(d["dim_sweep"], key=lambda r: r["samples_per_batch"])
    cmap = plt.get_cmap("viridis")

    fig, ax = plt.subplots(figsize=(7.8, 5.0))
    for i, r in enumerate(rows):
        h = r["history"]
        it = [x["iteration"] + 1 for x in h]
        en = [x["energy"] for x in h]
        hit = abs(en[-1] - exact) < 1e-7
        ax.plot(it, en, "o-", color=cmap(i / (len(rows) - 1)), lw=2.2,
                ms=7 if not hit else 8.5, zorder=4 + hit,
                label=f"{r['samples_per_batch']} samples/batch → "
                      f"dim {r['subspace_dim']:,}" + ("  ✓ exact" if hit else ""))
    ax.axhline(exact, color=C2, lw=2.2, ls="--", zorder=6)
    ax.text(0.95, exact, f"exact  {exact:.6f} Ha", color=C2, fontsize=9.5,
            va="top", ha="left", fontweight="semibold")

    ax.set_xlabel("configuration-recovery iteration")
    ax.set_ylabel("energy (Hartree)")
    ax.set_title("Energy vs configuration-recovery iteration\n"
                 "(10e,10o), 20 qubits, 63,504 determinants")
    ax.set_xticks(range(1, 11))
    ax.set_xlim(0.7, 10.4)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.2f}"))
    ax.legend(loc="upper right", fontsize=8.6)
    fig.tight_layout()
    p = FIGS / "energy_vs_iteration_20q.png"
    fig.savefig(p)
    plt.close(fig)
    print(f"  wrote {p.name}")


if __name__ == "__main__":
    load()
    fig_runtime_vs_active_space()
    fig_gpu_sectors()
    fig_energy_vs_iteration_20q()
    print("done")

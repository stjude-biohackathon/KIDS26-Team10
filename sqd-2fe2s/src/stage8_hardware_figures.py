"""Stage 8 - figures for the real-quantum-hardware runs.

Everything here is plotted from the JSON that `stage7_hardware.py` wrote, so no QPU
time and no re-simulation is needed to regenerate the panels:

    python3 src/stage8_hardware_figures.py

Three figures:

  hw_yield_vs_depth.png   usable-shot fraction vs transpiled two-qubit gate count
  hw_energy_by_size.png   energy error after post-processing, per active space
  hw_mitigation_ab.png    usable-shot fraction with and without DD + twirling
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


def load():
    """Collect every hardware run into one list of dicts.

    The 20-qubit run predates the --tag/--suffix flags, so its result file was later
    overwritten by the 32-qubit run.  Its numbers are recovered from the consolidated
    `hardware_all.json` instead; the raw samples for it are still on disk
    (`hardware_samples_sto-3g_fe3d.json`) if anyone wants to re-derive them.
    """
    allj = json.load(open(RESULTS / "hardware_all.json"))
    runs = []
    for r in allj["runs"]:
        runs.append(dict(
            qubits=r["qubits"], space=r["active_space"], job=r["job_id"],
            gates=r["two_qubit_gates"], yield_n=r["frac_valid_N"],
            yield_nsz=r["frac_valid_N_and_Sz"], energy=r.get("energy"),
            reference=r.get("reference"), err=r.get("error_mha"),
            mitigated=bool(r.get("mitigation")),
        ))
    return allj, runs


def fig_yield_vs_depth(runs):
    """Usable-shot yield against circuit size - the result that actually generalises.

    'Usable' means the shot came back with both the right electron count and the right
    spin projection; those are the only shots the classical solver can use.  Plotted
    on a log axis because the fall is two orders of magnitude.
    """
    plain = sorted([r for r in runs if not r["mitigated"]], key=lambda r: r["qubits"])
    g = np.array([r["gates"] for r in plain], float)
    yn = np.array([r["yield_n"] for r in plain]) * 100
    ys = np.array([r["yield_nsz"] for r in plain]) * 100
    q = [r["qubits"] for r in plain]

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.plot(g, yn, "o--", color=INK3, lw=1.8, ms=8, label="correct electron count")
    ax.plot(g, ys, "o-", color=C1, lw=2.8, ms=10, zorder=5,
            label="correct count AND spin  (the shots SQD can use)")
    ax.set_yscale("log")
    for gi, yi, qi in zip(g, ys, q):
        ax.annotate(f"{qi} qubits\n{yi:.2f}% usable", (gi, yi),
                    textcoords="offset points", xytext=(14, -6), ha="left",
                    fontsize=9.5, color=C1, fontweight="semibold")
    ax.text(0.30, 0.13,
            f"{ys[0] / ys[-1]:.0f}x fewer usable shots,  {g[-1] / g[0]:.1f}x the gates",
            transform=ax.transAxes, color=C2, fontsize=10, ha="center")
    ax.set_xlabel("two-qubit gates in the transpiled circuit")
    ax.set_ylabel("fraction of 100,000 shots (%, log scale)")
    ax.set_title("Usable-shot fraction vs transpiled circuit size\n"
                 "IBM ibm_fez (156-qubit Heron), 100,000 shots per point")
    ax.set_ylim(0.05, 60)
    ax.set_xlim(1700, 7000)
    ax.legend(loc="upper right")
    fig.tight_layout()
    p = FIGS / "hw_yield_vs_depth.png"
    fig.savefig(p)
    plt.close(fig)
    print(f"  wrote {p.name}")


def fig_energy_by_size(runs):
    """Energy error after post-processing, hardware, by problem size.

    Only the 20-qubit point has an exact reference that is genuinely exact; the
    36-qubit reference is extrapolated from the spin ladder and is marked as such.
    """
    plain = sorted([r for r in runs if not r["mitigated"]], key=lambda r: r["qubits"])
    labels = [f"{r['qubits']} qubits\n{r['space']}" for r in plain]
    errs = [r["err"] if r["err"] is not None else np.nan for r in plain]
    x = np.arange(len(plain))

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    cols = [C3, C2, C2]
    bars = ax.bar(x, np.maximum(errs, 1e-3), color=cols, width=.58, zorder=3)
    ax.axhline(1.6, color=C7, lw=1.8, ls="--", zorder=4)
    # left edge: on the right it sat on top of the 36-qubit bar
    ax.text(-0.30, 1.9, "chemical accuracy 1.6 mHa", color=C7, fontsize=9, ha="left")
    ax.set_yscale("log")
    for xi, e, b in zip(x, errs, bars):
        txt = "exact\n(+0.0000 mHa)" if e < 1e-6 else f"+{e:.0f} mHa"
        ax.annotate(txt, (xi, max(e, 1e-3)), textcoords="offset points",
                    xytext=(0, 6), ha="center", fontsize=9.5, color=INK,
                    fontweight="semibold")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("energy error vs exact diagonalisation (mHa, log scale)")
    ax.set_ylim(5e-4, 4000)
    ax.set_title("Hardware energy error after classical post-processing\n"
                 "IBM ibm_fez, 100,000 shots per active space")
    ax.text(.02, .78, "36-qubit reference is extrapolated\nfrom the exact spin ladder",
            transform=ax.transAxes, fontsize=8.5, color=INK2)
    fig.tight_layout()
    p = FIGS / "hw_energy_by_size.png"
    fig.savefig(p)
    plt.close(fig)
    print(f"  wrote {p.name}")


def fig_mitigation_ab(runs):
    """Dynamical decoupling + twirling, on and off, at 32 and 36 qubits.

    The metric that matters is the usable fraction (count AND spin).  It does not move.
    """
    # only the sizes where BOTH arms actually ran
    have = {(r["qubits"], r["mitigated"]) for r in runs}
    groups = [q for q in sorted({r["qubits"] for r in runs})
              if (q, True) in have and (q, False) in have]
    if not groups:
        print("  skipped hw_mitigation_ab.png - no size has both arms on disk")
        return
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    width = 0.34
    x = np.arange(len(groups))
    for k, (key, col, lab) in enumerate([
            ("plain", INK3, "plain circuit"),
            ("mit", C1, "dynamical decoupling (XpXm) + gate/measure twirling")]):
        vals, ann = [], []
        for qb in groups:
            r = [z for z in runs if z["qubits"] == qb
                 and z["mitigated"] == (key == "mit")][0]
            vals.append(r["yield_nsz"] * 100)
            ann.append(r["gates"])
        ax.bar(x + (k - .5) * width, vals, width, color=col, zorder=3, label=lab)
        for xi, v, gg in zip(x + (k - .5) * width, vals, ann):
            ax.annotate(f"{v:.2f}%\n{gg:,} 2q gates", (xi, v),
                        textcoords="offset points", xytext=(0, 5), ha="center",
                        fontsize=9, color=INK)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{q} qubits" for q in groups])
    ax.set_ylabel("usable shots: correct electron count AND spin (%)")
    ax.set_ylim(0, 0.72)
    ax.set_title("Usable-shot fraction with and without error suppression\n"
                 "IBM ibm_fez, 100,000 shots per bar")
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    p = FIGS / "hw_mitigation_ab.png"
    fig.savefig(p)
    plt.close(fig)
    print(f"  wrote {p.name}")


if __name__ == "__main__":
    allj, runs = load()
    print(f"{len(runs)} hardware runs on {allj['backend']}")
    fig_yield_vs_depth(runs)
    fig_energy_by_size(runs)
    fig_mitigation_ab(runs)
    print("done")

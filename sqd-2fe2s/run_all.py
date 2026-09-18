#!/usr/bin/env python3
"""Run the whole pipeline: classical reference -> SQD -> figures and RESULTS.md."""

from __future__ import annotations

import argparse
import subprocess
import sys
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent
SRC = ROOT / "src"


def run(script: str, *extra: str) -> None:
    cmd = [sys.executable, str(SRC / script), *extra]
    print(f"\n$ {' '.join(cmd)}\n", flush=True)
    subprocess.run(cmd, check=True, cwd=ROOT)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--basis", default="sto-3g")
    ap.add_argument("--avas-threshold", type=float, default=0.5)
    ap.add_argument("--active", default="fe3d")
    ap.add_argument("--shots", type=int, default=100_000)
    ap.add_argument("--quick", action="store_true",
                    help="small sweeps, for checking the pipeline runs end to end")
    ap.add_argument("--optimize-ansatz", action="store_true")
    args = ap.parse_args()

    # One tag, passed explicitly to every stage.  Letting stage 1 fall back to its own
    # default while telling later stages a different name means they silently read a
    # stale stage-1 file instead of the one just written.
    tag = f"{args.basis}_{args.active}"

    run("stage1_classical.py", "--basis", args.basis,
        "--avas-threshold", str(args.avas_threshold),
        "--active", args.active, "--tag", tag)

    stage2 = ["--tag", tag, "--shots", str(args.shots)]
    if args.quick:
        stage2 += ["--shots", "20000", "--depol-sweep", "0.0", "0.01",
                   "--spb-sweep", "50", "200", "--max-iterations", "4"]
    if args.optimize_ansatz:
        stage2 += ["--optimize-ansatz"]
    run("stage2_sqd.py", *stage2)

    if not args.quick:
        run("stage2b_orbital_opt.py", "--tag", tag, "--shots", str(args.shots),
            "--oo-iters", "12")

    stage3 = ["--tag", tag]
    if args.optimize_ansatz:
        stage3 += ["--suffix", "_opt"]
    run("stage3_report.py", *stage3)

    print("\ndone - see results/RESULTS_*.md and results/figures/")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Collect the K=10 batch array into the single number SQD would report.

SQD keeps the lowest energy across batches, so the deliverable is min(), not mean().
The spread across seeds is itself informative: it shows how much of SQD's answer is
luck of the draw in which configurations got sampled.
"""
import glob, json, statistics

rows = [json.load(open(f)) for f in sorted(glob.glob("results/batches/task*.json"))]
if not rows:
    raise SystemExit("no results yet in results/batches/")
rows.sort(key=lambda r: r["err_mha"])
print(f"{'task':>5} {'seed':>6} {'subspace':>11} {'% CAS':>7} {'err (mHa)':>11} {'min':>7}")
for r in rows:
    print(f"{r['task']:>5} {r['seed']:>6} {r['subspace_dim']:>11,} "
          f"{r['pct_cas']:>6.2f}% {r['err_mha']:>11.3f} {r['wall_s']/60:>7.1f}")
errs = [r["err_mha"] for r in rows]
print(f"\nSQD reports the best batch : {min(errs):+.3f} mHa   <- the deliverable")
print(f"worst batch                : {max(errs):+.3f} mHa")
if len(errs) > 1:
    print(f"spread across {len(errs)} seeds     : {max(errs)-min(errs):.3f} mHa "
          f"(sd {statistics.stdev(errs):.3f})")
print(f"\nbaseline to beat (K=3, 200k shots, 10 iters, 4.54 h): +34.478 mHa")

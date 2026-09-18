# Cluster jobs — times, and what to cut when the clock runs out

Three **independent** jobs. Submit them separately (or all at once — they don't depend
on each other). Independence is the point: a slow rung cannot eat another rung's time,
which is what would happen if they were chained in one script.

```bash
bash  preflight.sh                  # run this first, takes under a minute
sbatch run_rungA.sbatch             # rung A: (22e,16o), 32 qubits - start here
sbatch run_casci_26e18o.sbatch      # rung B: (26e,18o), 36 qubits - the expensive one
sbatch run_rungC.sbatch             # rung C: (34e,22o), 44 qubits - parallel with A
```

**Status, September 2026: rungs A and B are DONE** - see
[../docs/RESULTS_MASTER.md](../docs/RESULTS_MASTER.md). Rung A ran on a laptop
(0.96 h CASCI, 4.54 h SQD) and rung B on an A100 GPU (4 of 6 spin sectors). These
scripts are kept so the runs are reproducible on a real cluster, and so rung C - which
was never run - has somewhere to start.

## Times

Extrapolated from **measured** per-sector timings on 8 cores (`t ≈ 1.45e-6 · dim^1.22`,
fitted on three measured points), so treat them as ±50%, not gospel. More cores helps
sub-linearly; memory is the real constraint.

| rung | active space | qubits | determinants | CASCI ladder | SQD | wall requested | memory |
|---|---|---|---|---|---|---|---|
| **A** | (22e,16o) +bridging S 3p | 32 | 19,079,424 | **30–40 min** | 1–2 h | 3 h | 32 GB |
| **B** | (26e,18o) +Fe 4s | 36 | 73,410,624 | **2–2.5 h** | 30–60 min | 6 h | 96 GB |
| **C** | (34e,22o) +2 thiolates | 44 | 693,479,556 | *not possible* | 1–1.5 h | 3 h | 64 GB |

Measured rung-A sectors (the basis for the fit):

| Sz | determinants | time |
|---|---|---|
| 5 | 8,008 | 0.1 s |
| 4 | 183,040 | 3.9 s |
| 3 | 1,544,400 | 52.8 s |
| 2 | 6,406,400 | ~5 min *(extrap.)* |
| 1 | 14,574,560 | ~14 min *(extrap.)* |
| 0 | 19,079,424 | ~19 min *(extrap.)* |

**All three rungs in parallel: ~3 h wall. Sequentially: ~5–6 h.**

## Plan B — what to cut, in order

The presentation is already safe: **the (10e,10o) result is banked and needs no cluster
at all.** Everything here is upside. Cut from the bottom of this list up.

| # | if this happens | do this | you still have |
|---|---|---|---|
| 0 | *nothing runs at all* | present as-is | Full validated result: SQD exact vs CASCI, noise sweep, controls, scaling crossover. Meets every criterion in the brief. |
| 1 | **queue is slow / < 4 h left** | submit **rung A only** | The J improvement — the single most valuable addition. A alone is the headline. |
| 2 | **rung B overruns** | let it die on the wall clock | The ladder is **checkpointed per sector**, cheapest first. A killed job still leaves `stage1_*_ladder_partial.json`, and **3 clean spin states are enough to fit J**. So even a half-finished rung B yields a usable number. |
| 3 | **rung C's sampling is too slow** | kill it; drop the 44-qubit point | Rungs A and B, which have exact references and are the stronger science anyway. |
| 4 | **rung A's SQD sweeps overrun** | kill after `[A1]` | The exact CASCI ladder, hence J at (22e,16o). That is the deliverable; the SQD sweeps are confirmation. |
| 5 | **everything overruns** | `--spb-sweep 200` only, `--depol-sweep 0.0` | One SQD point per rung instead of a sweep. Cheap and still a valid comparison. |

### Built-in safety

- **Cheapest sector first.** The Sz=0 sector is 2,400× bigger than Sz=5. Running it
  first would mean learning nothing until the hardest part finished; reversed, you get
  usable partial results within a minute.
- **Per-sector checkpointing** to `results/stage1_<tag>_ladder_partial.json`, with a J
  fit on whatever is complete so far.
- **Independent jobs**, so one overrun doesn't cascade.
- **Rung C's stage 1 is already done and verified locally** (95 s), so that job cannot
  fail on setup.

### The one thing not to attempt

**Do not try to get J out of SQD.** Measured on (10e,10o) against the known exact
answer: per-sector SQD errors reach 55 mHa while the whole ladder spans 3.95 mHa, giving
J = −174.6 cm⁻¹ against a true 29.0 — wrong sign. No amount of cluster time fixes it;
J is a difference of near-degenerate states and SQD delivers absolute energies. Every J
in these jobs comes from exact CASCI. Evidence: `results/stage4_sto-3g_fe3d.json`.

## What success looks like

J should climb toward the published DMRG value as the active space grows:

| active space | J (cm⁻¹) |
|---|---|
| (10e,10o) Fe 3d — **done** | **28.8** |
| (22e,16o) +bridging S — rung A | *expected to rise* |
| (26e,18o) +Fe 4s — rung B | *expected to rise further* |
| DMRG (30e,20o), TZP-DKH — published | 236 |
| Experiment | 148 ± 16 |

If rung A's J does *not* rise, that is still a reportable result — it would mean the
minimal basis, not the active space, is the dominant error. Either outcome is a slide.

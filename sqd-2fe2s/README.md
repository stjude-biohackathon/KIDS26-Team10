# SQD for iron–sulfur electronic structure — [Fe₂S₂(SMe)₄]²⁻

**Sample-Based Quantum Diagonalization on a [2Fe–2S] cluster, validated against the exact
classical answer, and run on real IBM quantum hardware.**

BioHackathon project, September 2026. A direct port of the
[IBM SQD chemistry tutorial](https://quantum.cloud.ibm.com/docs/en/tutorials/sample-based-quantum-diagonalization)
with N₂ swapped for an iron–sulfur cluster — plus the exact reference, the control
experiments, and the negative results that a quantum-chemistry claim needs.

> **📄 New here? Read [REPORT.md](REPORT.md).** It is the self-contained project report:
> problem, method, every result, every caveat, and the questions we expect to be asked.

---

## Results in one table

| active space | qubits | configurations | CASCI (exact) | SQD | **J (cm⁻¹)** |
|---|---|---|---|---|---|
| (10e,10o) Fe 3d | 20 | 63,504 | −5013.64906670 | **+0.0000 mHa** ✅ | **28.8** |
| (22e,16o) + bridging S 3p | 32 | 19,079,424 | −5013.72724223 | **+34.5 mHa** | **46.7** |
| (26e,18o) + Fe 4s | 36 | 73,410,624 | 4/6 spin sectors (GPU) | +653 mHa *(hardware)* | **82.0** * |
| DMRG (30e,32o), better basis — *published* | — | — | — | — | 236 |
| **Experiment** | — | — | — | — | **148 ± 16** |

\* four-sector fit, singlet extrapolated.

**Real hardware** — IBM `ibm_fez` (156-qubit Heron), five jobs, 100,000 shots each:

| qubits | 2q gates | correct N | **usable (N and Sz)** | energy error |
|---|---|---|---|---|
| **20** | 2,031 | 17.22% | **5.90%** | **+0.0000 mHa — exact** |
| 32 | 4,878 | 1.40% | **0.47%** | +529.7 mHa |
| 36 | 5,982 | 0.31% | **0.10%** | +653.4 mHa |

At 20 qubits **94% of the shots came back physically impossible** and the classical
post-processing still recovered the exact energy in the correct spin state. That is the
result — with the qualification that the recovered subspace was 100% of the space, so it
demonstrates noise resilience, not compression.

### And the honest part

- On this molecule SQD needs **~30% of the whole configuration space** to reach 34.5 mHa.
  The wavefunction is not compact (participation ratio ~1,956).
- At 32 qubits **exact diagonalization was 4.7× faster *and* exact** (0.96 h vs 4.54 h).
- At 20 qubits a **uniform-random-bitstring control matches the quantum circuit**.
- **J cannot be extracted from SQD at all** — it comes out with the wrong sign.

All four are measured, not asserted. An independent 2026 assessment of SQD
([arXiv:2501.07231](https://arxiv.org/abs/2501.07231)) reaches the same conclusion on the
same systems.

---

## Quick start

```bash
pip install -r requirements.txt
python3 run_all.py
```

~12 minutes on a laptop, and reproduces the entire validated (10e,10o) result: classical
reference, exact CASCI spin ladder, SQD with noise and subspace sweeps, both control
experiments, figures and a markdown report into `results/`.

Smoke test in ~2 minutes:

```bash
python3 run_all.py --quick
```

---

## The pipeline

```
src/geometry.py           idealised [Fe2S2(SMe)4]2-  (24 atoms, charge -2)
  |
src/stage1_classical.py   RHF (DF + newton) -> AVAS -> active-space integrals
                          -> exact CASCI spin ladder S = 0..5  (ground truth + J)
                          -> frozen-core CCSD -> t1/t2 amplitudes
  |
src/stage2_sqd.py         LUCJ circuit -> sample -> post-select
                          -> self-consistent configuration recovery
                          -> subspace diagonalisation
                          + noise sweep, subspace sweep, and two controls
  |
src/stage3_report.py      figures + results/RESULTS_<tag>.md
```

Side branches:

| script | what it does |
|---|---|
| `src/stage2b_orbital_opt.py` | orbital optimisation at frozen subspace — recovers only 3–16% of the gap |
| `src/stage4_spin_ladder.py` | SQD per spin sector + Heisenberg fit. **A documented failure** — wrong sign for J. Kept because the negative result is the finding |
| `src/stage5_scaling.py` | nested-slice timing study, 6 → 16 orbitals: where CASCI hits its wall |
| `src/stage6_talk_figures.py` | the three presentation figures |
| `src/stage7_hardware.py` | **real QPU runs** on IBM hardware — see [docs/HARDWARE.md](docs/HARDWARE.md) |
| `src/collect_hardware.py` | consolidates every QPU run into `results/hardware_all.json` |
| `src/stage8_hardware_figures.py` | the three hardware figures |
| `src/stage9_extra_figures.py` | coverage curve, CASCI-vs-SQD cost, and the shot-count null result |
| `src/stage10_timing_figures.py` | **runtime vs active space** (measured), GPU sector scaling, and energy-in-Hartree vs iteration |
| `src/ffsim_patch.py` | **works around an upstream `ffsim` bug** that silently breaks the tutorial's `optimize=True`. See below |

Individual stages:

```bash
python3 src/stage1_classical.py --basis sto-3g --avas-threshold 0.5 --active fe3d
python3 src/stage2_sqd.py       --tag sto-3g_fe3d --shots 100000
python3 src/stage3_report.py    --tag sto-3g_fe3d
```

Active spaces available to `--active`: `fe3d`, `fe3d+brs3p`, `fe3d+brs3p+fe4s`,
`fe3d+brs3p+2terms3p`, `fe3d+alls3p`.

---

## Real quantum hardware

Five completed jobs on `ibm_fez`, ~250 s of a 600 s allocation. Full procedure, including
the credential setup and the mistakes that cost budget, is in
**[docs/HARDWARE.md](docs/HARDWARE.md)**.

```bash
python3 src/stage7_hardware.py --tag sto-3g_fe3d --shots 100000 --estimate  # free
python3 src/stage7_hardware.py --tag sto-3g_fe3d --shots 100000             # submits
```

**Nothing in this repository contains or reads an API token.** `stage7_hardware.py` calls
`QiskitRuntimeService()` with no arguments and picks up the account you saved yourself with
`QiskitRuntimeService.save_account(...)` — see docs/HARDWARE.md §1.

**Every job's raw shots are committed** (`results/hardware_samples_*.json`), so all
post-processing replays at zero QPU cost:

```bash
python3 src/stage7_hardware.py --tag sto-3g_fe3d --replay
```

---

## Simulator backends

No QPU needed for anything except stage 7. `stage2_sqd.py --backend`:

| backend | what it is | cost (20 qubits, 3k shots) |
|---|---|---|
| `ffsim` *(default)* | exact statevector sampling inside the fixed particle-number sector, with an optional global depolarizing channel | ~1 s |
| `aer` | `qiskit-aer` statevector with a gate-level depolarizing noise model (1-qubit `depol`, 2-qubit 10×) — the closest stand-in for a device | ~13 s |
| `aer-mps` | the same, matrix-product-state method; for when the qubit count outgrows a dense statevector | — |

`ffsim` is the default because it is ~10× faster *and* because it samples only physical
(correct-N) states, so the depolarizing knob is the only source of invalid shots and the
sampling-efficiency numbers stay interpretable.

⚠️ **Aer must transpile at `optimization_level=0`.** Levels ≥ 1 hang on ffsim's custom LUCJ
gate.

---

## Measured runtimes

Laptop, 8 cores, 16 GB, STO-3G:

| what | time |
|---|---|
| **`run_all.py` — full (10e,10o) pipeline** | **~12 min** |
| `run_all.py --quick` | ~2 min |
| (22e,16o) CASCI — 19.1M determinants, 1 core, idle machine | **0.96 h** |
| (22e,16o) SQD — 10 recovery iterations → 5.7M subspace | **4.54 h** |
| (26e,18o) CASCI, 4 spin sectors | several h on an **A100 GPU** (needs 17.6 GB) |
| one hardware job, 100k shots | ~54 s of QPU time |

Both 32-qubit timings were measured on an **otherwise idle machine**. Earlier figures were
contaminated by concurrent jobs badly enough to reverse the conclusion — if you re-measure,
measure alone.

Larger active spaces on a cluster: run `bash cluster/preflight.sh` first (it takes under a
minute and catches the OS/BLAS/Python differences cheaply), then see
[cluster/README.md](cluster/README.md) for the rungs and the cut list.

---

## Three things to know before reading the results

**1. An upstream `ffsim` bug blocks the tutorial's `optimize=True`.**
`ffsim.linalg.util.unitaries_to_parameters` calls `scipy.linalg.logm(mats)` on a batch of
shape `(n, dim, dim)`, but `logm` accepts only one square matrix →
`ValueError: expected square array_like input`. That breaks `UCJOp*.to_parameters()` and
therefore `from_t_amplitudes(..., optimize=True)`, which the official tutorial uses.
Reproduced on ffsim 0.0.83 **and** 0.0.84. `src/ffsim_patch.py` maps the logarithm over the
batch; the round-trip is then exact (`|⟨ψ|ψ′⟩| = 1.000000000000`). It installs only if the
bug is present, so it no-ops once upstream fixes it. Measured impact: **+329 mHa → +0.0000
mHa**.

**2. A silent PySCF Davidson trap.** In the near-degenerate Sz = 0 sector `max_space=12`
(the default) stops short of convergence — and asking for *more roots makes it worse*
(+0.022 mHa at `nroots=1`, +0.128 mHa at `nroots=3`). `max_space=30` is exact and faster.
The S=0→1 gap is 0.245 mHa, so the default biases J by ~50%. `lowest_in_sector()` sets it.

**3. The (10e,10o) space is exact but small.** 63,504 determinants means CASCI is cheap and
SQD has almost nothing to compress. Two consequences are reported rather than hidden: the
random-bitstring control performs as well as the circuit, and the whole spin ladder spans
under 4 mHa, which makes the usual 1.6 mHa chemical-accuracy target uninformative here.

---

## Classical baselines

Implemented and run:

- **CASCI / FCI in the active space** — at these sizes this *is* full CI, i.e. the exact
  solution of the same Hamiltonian SQD attacks. The ground truth every number is measured
  against.
- **Selected CI** — diagonalise in the product closure of the largest-|c| CASCI
  determinants. The closest classical analogue of SQD. A comparison strategy, *not* an
  upper bound: at small subspace dimensions the circuit's sampled subspace actually beats
  it.
- **RHF and frozen-core CCSD** — the mean-field and single-reference baselines. CCSD **does
  not converge** on this half-filled d shell, which is itself the diagnostic that the
  system is outside single-reference range.
- **Heisenberg fit** to the exact CASCI S = 0…5 ladder → the exchange coupling J, the
  quantity comparable to experiment. Convention `H = 2 J S₁·S₂`, `E(S) = J S(S+1)`,
  `J > 0` antiferromagnetic — the same as the reference literature.

Not implemented (stretch goals; only *meaningful* once the active space outgrows CASCI):
DMRG via `block2`, SHCI via Dice, FCIQMC. `qiskit-addon-sqd` ships a Dice integration
guide if you want SHCI as a scalable cross-check.

---

## Repository layout

```
REPORT.md                   the project report - start here
README.md                   this file
requirements.txt
run_all.py                  one-command driver

src/
  geometry.py               cluster builder + bond-length sanity report
  common.py                 paths, constants, JSON helpers
  ffsim_patch.py            upstream-bug workaround
  stage1_classical.py       classical reference and exact ground truth
  stage2_sqd.py             the SQD runs (simulator)
  stage2b_orbital_opt.py    orbital optimisation at frozen subspace
  stage3_report.py          figures + RESULTS_<tag>.md
  stage4_spin_ladder.py     J from SQD - a documented failure
  stage5_scaling.py         CASCI vs SQD timing study
  stage6_talk_figures.py    presentation figures
  stage7_hardware.py        real QPU runs
  collect_hardware.py       consolidate QPU results
  stage8_hardware_figures.py  hardware figures
  stage9_extra_figures.py   coverage / cost / shot-count figures
  stage10_timing_figures.py runtime-vs-active-space and energy-vs-iteration figures

docs/
  RESULTS_MASTER.md         single source of truth for every number
  HARDWARE.md               how the QPU runs were done, step by step
  SUBMISSION.md             pitch, abstract, methods
  PROPOSAL.md               narrative version, no orbital jargon
  EXPLAINER.md              lay-audience explanation
  PLAN.md                   the agreed plan, task split, pitfalls
  RUN_SUMMARY.md            what was run, with timings
  LITERATURE_COMPARISON.md  per-active-space placement vs published work
  REFERENCES.md             every paper and package, and what it was for

cluster/                    sbatch scripts + preflight check + cut list
notebook/                   GPU notebook for the (22e,16o) space
geometry/                   generated .xyz structure
results/                    JSON, raw QPU shots, reports, figures/
```

---

## Requirements

```bash
pip install -r requirements.txt
```

Tested on Python 3.13 with pyscf 2.14.0, qiskit 2.4.2, qiskit-addon-sqd 0.13.1,
ffsim 0.0.84. Hardware runs additionally need `qiskit-ibm-runtime`; the GPU CASCI path
needs `gpu4pyscf`.

---

## References

Full list with what each was used for: [docs/REFERENCES.md](docs/REFERENCES.md). The two
that matter most:

- **Sharma, Sivalingam, Neese & Chan**, *Nature Chemistry* **6**, 927 (2014),
  [arXiv:1408.5080](https://arxiv.org/abs/1408.5080) — the classical benchmark, same model
  complex, DMRG J = 236 cm⁻¹.
- **Robledo-Moreno et al.**, [arXiv:2405.05068](https://arxiv.org/abs/2405.05068) — SQD on
  the same two clusters at 45 and 77 qubits, on 3,072 Fugaku cores. The quantum benchmark.

---

## License

MIT — see [LICENSE.md](../LICENSE.md) at the repository root.

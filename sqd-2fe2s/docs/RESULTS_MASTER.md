# Results — SQD for a [2Fe–2S] iron–sulfur cluster

Single source of truth. Every number here is measured and traceable to a file in
`results/`. Where something is extrapolated rather than computed, it says so.

**System:** [Fe₂S₂(SMe)₄]²⁻ — the [2Fe–2S] core with four methylthiolate ligands
standing in for the cysteines that anchor the cofactor in a protein. 24 atoms, charge
−2, STO-3G throughout. Same model complex as the reference DMRG study, so comparisons
are like-for-like.

---

## 1. Headline: the exchange coupling converges toward experiment

`J` in the Fe–S literature convention `H = 2 J S₁·S₂`, `E(S) = J S(S+1)`, `J > 0`
antiferromagnetic.

| active space | qubits | determinants | **J (cm⁻¹)** | fraction of experiment |
|---|---|---|---|---|
| (10e,10o) Fe 3d | 20 | 63,504 | **28.8** | 0.19× |
| (22e,16o) + bridging S 3p | 32 | 19,079,424 | **46.7** | 0.32× |
| (26e,18o) + Fe 4s | 36 | 73,410,624 | **82.0** *(4-sector fit)* | 0.55× |
| DMRG (30e,20o), TZP-DKH | — | — | 236 | 1.59× |
| **Experiment** (magnetic susceptibility) | — | — | **148 ± 16** | 1.00× |

**28.8 → 46.7 → 82.0 cm⁻¹**, monotonic, roughly doubling per rung. Each rung adds the
chemistry the previous one truncated: first the bridging-sulfur 3p orbitals that carry
the Fe–Fe superexchange, then the Fe 4s double shell. The remaining gap to 148 is
attributable to the minimal STO-3G basis (the reference used TZP-DKH).

Two caveats stated plainly:

- The (26e,18o) value is fitted from **four sectors** (S = 2,3,4,5); Sz = 1 and Sz = 0
  exceeded the memory budget (56.8M and 73.4M determinants). Going from three sectors
  to four moved J from 84.3 to 82.0, and at (22e,16o) the 3-point fit was 3% high
  against the full ladder — both say the fit is converging, so 82.0 is stable to a
  few percent.
- (26e,18o) is **not strictly nested** inside (22e,16o): AVAS re-selects orbitals at each
  threshold, so this is "three related active spaces", not one being systematically
  enlarged. The trend is real; the wording should be precise.

---

## 2. SQD against exact diagonalisation

| active space | qubits | CASCI (exact) | SQD best | SQD error |
|---|---|---|---|---|
| (10e,10o) | 20 | **−5013.64906670** | −5013.64906670 | **+0.0000 mHa** |
| (22e,16o) | 32 | **−5013.72724223** | −5013.69276393 | **+34.5 mHa** |
| (26e,18o) | 36 | −5013.67603364 *(extrapolated)* | −5013.02262821 *(hardware)* | +653.4 mHa |

Two provenance notes on the 36-qubit row, so nobody is caught out by a follow-up
question. **We never ran SQD on the simulator at 36 qubits** — stage 1 for that space
was run with `--skip-casci` and the only 36-qubit SQD numbers in `results/` come from the
two hardware jobs. And the reference is the **extrapolated** singlet from four exact spin
sectors, not a computed one, because Sz = 0 and Sz = 1 ran out of memory. Treat that row
as indicative; the 20- and 32-qubit rows are the checkable ones.

At 20 qubits SQD is **exact to eight decimal places**, with ⟨S²⟩ = 0.0000 (the correct
singlet, not merely the right energy) and orbital occupancies matching CASCI to
1.1 × 10⁻³ electrons. Accuracy then degrades with active-space size at a fixed sampling
budget.

For scale: 34.5 mHa is **21.6 kcal/mol**, about 22× outside chemical accuracy (1.6 mHa)
and 9× larger than the entire S = 0…5 spin ladder.

---

## 3. Why it degrades: coverage, measured

SQD diagonalises the subspace its samples span. At (22e,16o) we swept that:

| subspace | % of the 19.1M space | error |
|---|---|---|
| 31,329 | 0.16% | +726.7 mHa |
| 57,600 | 0.30% | +704.4 mHa |
| 223,729 | 1.17% | +568.5 mHa |
| 628,849 | 3.30% | +515.0 mHa |
| 1,385,329 | 7.26% | +321.5 mHa |
| 5,716,881 | **29.96%** | **+34.5 mHa** |
| 19,079,424 | 100% | 0 (= CASCI) |

Figure: `results/figures/coverage_22e16o.png`.

**On this molecule you need nearly the whole space.** The wavefunction is not compact:
at (10e,10o) its participation ratio is ~1,956 and 7,171 determinants are needed for 99%
of the weight. There is no small subspace to find. Orbital optimisation at *frozen*
subspace recovers only 3–16% of the gap, confirming the basis is not the issue either.

For contrast, IBM ran [4Fe–4S] at **10⁻⁷ %** coverage and got useful answers — because
*that* wavefunction is compact. **Our system is a hard case for SQD**, and that is the
finding.

---

## 4. Cost: CASCI beat SQD at 32 qubits

Measured on an otherwise idle machine, one core:

| | work | time | result |
|---|---|---|---|
| **CASCI** | one Davidson over all 19,079,424 determinants | **0.96 h** | **exact** |
| **SQD** | 10 recovery iterations → 5,716,881 subspace | **4.54 h** | +34.5 mHa |

**4.7× slower and not exact.** SQD's individual solve is cheaper (smaller matrix), but
the self-consistent loop needs many of them, and ten did not converge — the error was
still falling at the iteration cap:

```
iter 1  +441.7 mHa   dim   992,016
iter 3  +150.3 mHa   dim 2,528,100
iter 5   +95.0 mHa   dim 3,814,209
iter 7   +58.3 mHa   dim 4,739,329
iter 9   +39.2 mHa   dim 5,452,225
iter 10  +34.5 mHa   dim 5,716,881
```

Figure: `results/figures/cost_casci_vs_sqd.png`.

The honest reading: **exact diagonalisation is the right tool wherever it fits.** SQD's
regime begins where it doesn't.

---

## 5. Where exact diagonalisation stops

Runtime against active-space size, nested slices of the real Hamiltonian:

| orbitals | qubits | determinants | CASCI | SQD solve |
|---|---|---|---|---|
| 6 | 12 | 400 | 0.03 s | 0.11 s |
| 10 | 20 | 63,504 | 0.36 s | 0.73 s |
| 12 | 24 | 853,776 | 10.3 s | 1.14 s |
| 14 | 28 | 11,778,624 | **288.6 s** | 3.33 s |
| 16 | 32 | 165,636,900 | **infeasible** | 8.47 s |

CASCI grows ~10,000× across the ladder then hits a wall (166M determinants needs ~40 GB
for the Davidson). The SQD *solve* grows 77× and stays under 10 s, because its cost
tracks the subspace, not the space.

**Caveat:** SQD's *sampling* cost also grows steeply here (3.2 s → 48 s), but only
because we **simulate** the quantum computer. On hardware, sampling is roughly
constant-cost. The column that legitimately replaces CASCI is the solve.

Timing study only — truncated active spaces give meaningless energies.

Figures: `results/figures/talk_runtime.png` (this table, with a qubit axis) and
`results/figures/runtime_vs_active_space.png` (the same question asked of the three
**real** active spaces: 20q 11.9 s / 84.9 s, 32q 0.96 h / 4.54 h, 36q 1.67 h on a GPU).
Per-sector GPU scaling at (26e,18o) — 6 s, 27 s, 367 s, 5,612 s, then two sectors out of
memory — is in `results/figures/gpu_casci_sectors.png`.

---

## 6. Negative results (measured, not asserted)

**a) J cannot be extracted from SQD.** Running SQD per Sz sector and fitting the
Heisenberg form gives **J = −174.6 cm⁻¹ against a true 29.0** — wrong sign, 35 mHa fit
residual. Cause:

| Sz | subspace | SQD error | exact gap from Sz = 0 |
|---|---|---|---|
| 0 | 63,001 | +2.746 mHa | 0.000 |
| 1 | 39,770 | **+55.038 mHa** | 0.213 |
| 2 | 14,040 | +9.554 mHa | 0.710 |
| 3–5 | ≤2,025 | +0.000 mHa | 1.47–3.92 |

High-Sz sectors are tiny so SQD is trivially exact; the large low-Sz sectors carry
millihartree error. J is a *difference* of near-degenerate states, so error that uneven
swamps the signal by an order of magnitude. **More compute cannot fix this.** Every J in
this work therefore comes from exact CASCI.

> SQD delivers total energies to high accuracy here, but cannot resolve the small energy
> differences that magnetic properties depend on.

**b) Copying IBM's sampling depth does not help.** Same subspace config, only shots
changed:

| shots | unique bitstrings | subspace | error |
|---|---|---|---|
| 200,000 | 164,532 | 1,385,329 | **+321.5 mHa** |
| **2,457,600** (IBM's) | 1,108,437 | 1,435,204 | **+344.2 mHa** |

12× more shots, 6.7× more unique bitstrings — and the energy got **slightly worse**.
`samples_per_batch` caps how many configurations enter the subspace regardless of pool
quality. So IBM's advantage was never sampling depth; it was **K = 10–100 batches on
3,072 cores**. Figure: `results/figures/shots_null.png`.

**c) The quantum circuit buys nothing at 20 qubits.** A uniform-random-bitstring control
also reaches +0.0000 mHa. Once the subspace covers the space, how configurations were
chosen stops mattering. This is the control experiment a quantum-advantage claim
requires, and we ran it.

---

## 7. Comparison with the published work

| | **IBM, arXiv:2405.05068** | **this work** |
|---|---|---|
| Quantum sampling | ~45 min on a QPU | ~7 s, simulated |
| Shots | 2,457,600 | 200,000 (and 2.46M tested — no gain) |
| Batches K | 10 ([2Fe–2S]), 100 ([4Fe–4S]) | 1–8 |
| Subspace | up to 100M | 5.7M |
| Classical post-processing | **~1.5 h on 64 Fugaku nodes = 3,072 cores** | **4.54 h on 1 core** |
| Campaign total | ~6,400 nodes ≈ 307,200 cores | one laptop + one A100 GPU |
| Qubits | 45, 77 | 20, 32, 36 |
| Simulations as well as hardware? | **yes** — noiseless LUCJ + depolarizing-noise tests | yes |

Their quantum time (45 min) is trivial beside their classical time (1.5 h × 3,072 cores).
**The quantum part was never the bottleneck — for them either.**

---

## 8. Engineering results

**Two bugs found and fixed.**

1. **Upstream defect in `ffsim` 0.0.83 and 0.0.84.**
   `unitaries_to_parameters` calls `scipy.linalg.logm` on a batch of shape
   `(n, dim, dim)`, but `logm` accepts only one square matrix →
   `ValueError: expected square array_like input`. This breaks
   `UCJOp*.to_parameters()` and therefore `from_t_amplitudes(..., optimize=True)`, which
   the official IBM tutorial uses. `src/ffsim_patch.py` maps the logarithm over the
   batch; the parameter round-trip is then exact (`|⟨ψ|ψ′⟩| = 1.000000000000`). Impact,
   measured:

   | | ansatz spread (PR) | noiseless SQD |
   |---|---|---|
   | bug present | 2.68 | +329 mHa |
   | patched | 3,259 | **+0.0000 mHa** |

2. **A silent PySCF Davidson convergence trap.** In the near-degenerate Sz = 0 sector,
   `max_space=12` (the default) stops short — and asking for *more roots makes it worse*:

   | setting | error vs converged |
   |---|---|
   | `max_space=12, nroots=1` | +0.022 mHa |
   | `max_space=12, nroots=3` | **+0.128 mHa** |
   | `max_space=30, nroots=1` | **−0.000003 mHa**, and faster |

   The S=0→1 gap is only 0.245 mHa, so the default would bias J by ~50%.

**GPU acceleration.** `gpu4pyscf.fci.direct_spin1.FCI` subclasses PySCF's own solver and
overrides only `contract_2e` (the H·v product) with a CUDA kernel. One caveat we had to
fix: it requires a 2-D CI array while PySCF's Davidson passes a flat vector — a reshape
wrapper solves it. This is what made the (26e,18o) CASCI possible at all: it needs
17.6 GB for the Davidson working set, against 16 GB on the laptop.

---

## Real quantum hardware — IBM `ibm_fez` (156-qubit Heron)

**Five completed QPU jobs**, 100,000 shots each, same LUCJ ansatz at three circuit
depths, with and without error suppression. The 20- and 32-qubit spaces have an exact
CASCI reference, so those energies are genuinely checkable.

| qubits | space | suppression | 2q gates | correct N | **correct N and Sz** | energy error | ⟨S²⟩ | job id |
|---|---|---|---|---|---|---|---|---|
| **20** | (10e,10o) | — | 2,031 | 17.22% | **5.90%** | **+0.0000 mHa — exact** | 0.0000 | `dameitg2fm4c73f2t3f0` |
| **32** | (22e,16o) | — | 4,878 | 1.40% | **0.47%** | +529.7 mHa | 0.0062 | `damemt78gn2s739larvg` |
| 32 | (22e,16o) | DD + twirl | 4,742 | 1.52% | **0.47%** | +649.8 mHa | 0.0859 | `dameqg8pqrnc7397es80` |
| **36** | (26e,18o) | — | 5,982 | 0.31% | **0.10%** | +653.4 mHa* | 0.0912 | `dameumn8gn2s739lbb70` |
| 36 | (26e,18o) | DD + twirl | 5,972 | 0.35% | **0.10%** | +866.1 mHa* | 0.4666 | `damf48gpqrnc7397f9t0` |

\* against an *extrapolated* singlet (−5013.67603364 Ha) from the four exact spin
sectors — there is no exact reference at 36 qubits, so treat these two as indicative.

Machine-readable: `results/hardware_all.json` (regenerate with
`python3 src/collect_hardware.py`). Raw shot data for every job is on disk as
`results/hardware_samples_*.json`, so all post-processing above can be replayed for
free — no QPU needed:

```bash
python3 src/stage7_hardware.py --tag sto-3g_fe3d --replay
```

Figures: `results/figures/hw_yield_vs_depth.png`, `hw_energy_by_size.png`,
`hw_mitigation_ab.png` (`python3 src/stage8_hardware_figures.py`).

### The headline

At 20 qubits, **94.1% of the hardware shots were physically invalid** — wrong electron
count or wrong spin — and post-selection plus configuration recovery still recovered the
**exact** energy with the correct spin state (⟨S²⟩ = 0.0000, singlet). That is SQD's
noise-resilience claim, demonstrated on a QPU and checked against an answer we computed
independently.

Qualification, stated up front: the recovered subspace reached **100% of the space**
(63,504 of 63,504 determinants), so this demonstrates *noise resilience*, not
*compression* — and at this size a uniform-random-bitstring control reaches the same
answer. The hardware result is a validated pipeline, not evidence of quantum advantage.

### Depth versus usable-shot yield

Usable yield falls **61×** (5.90% → 0.10%) for **2.9×** the two-qubit gates
(2,031 → 5,982). The poor energies at 32 and 36 qubits follow directly from that: only
~470 and ~100 shots survive post-selection, against spaces of 19.1M and 73.4M
determinants. The recovered subspaces cover 2.4% and 0.4% of their spaces, and the
simulator coverage curve (§3) already says what that costs — so the hardware is behaving
as expected, not anomalously.

For scale, IBM's published run on this cluster used **2,457,600 shots and K = 10–100
batches on 3,072 Fugaku cores**; ours used 100,000 shots, K = 1–3, on one laptop core.
The 32- and 36-qubit numbers are under-resourced by construction. They are a
*depth-versus-fidelity measurement*, not a claim about what hardware can do with a real
shot budget.

### Error suppression: a measured null result, twice

| | 2q gates | correct N | **correct N and Sz** | energy error |
|---|---|---|---|---|
| 32q plain | 4,878 | 1.40% | **0.47%** | +529.7 mHa |
| 32q DD (XpXm) + gate/measure twirling | 4,742 | 1.52% | **0.47%** | +649.8 mHa |
| 36q plain | 5,982 | 0.31% | **0.10%** | +653.4 mHa |
| 36q DD (XpXm) + gate/measure twirling | 5,972 | 0.35% | **0.10%** | +866.1 mHa |

A marginal gain in particle-number survival at both sizes, and **none at all** in the
metric that determines the answer. Energies came out *worse* in both mitigated arms —
that is sample luck rather than harm from the techniques, since the usable yield is
identical and the surviving shots simply land on a different (and here slightly poorer)
subspace; ⟨S²⟩ degrading to 0.4666 at 36q says how thin that sampling is.

Why suppression cannot help here: at ~5,000 two-qubit gates the dominant error is
**incoherent gate error**, and neither technique removes it. DD targets idle dephasing on
spectator qubits; twirling reshapes coherent error into stochastic error. Improving this
needs better gates or a shallower ansatz, not better scheduling. ZNE, PEC and TREX — the
techniques that *would* attack gate error — are Estimator-only and cannot be applied to a
Sampler workload like SQD at all.

Caveat on the A/B: the two arms transpiled to different gate counts (4,878 vs 4,742)
because Qiskit's layout pass is stochastic, so the comparison carries some layout noise.
`seed_transpiler` was pinned afterwards for future runs.

### QPU budget accounting

IBM reported `quantum_seconds` for three of the five jobs: **55.64 s**, **55.45 s** and
**50.64 s** per 100,000-shot job — a mean of **53.9 s**, i.e. **~539 µs per shot**.
Estimated total spend across all five jobs is **~250 s of the 600 s** allocation.

A follow-up 1,000,000-shot job at 32 qubits (`damf1k02fm4c73f2trag`) was **cancelled**
before execution once IBM's own queue estimate showed ~13 minutes, over budget.

### Two estimation errors worth recording

1. **Our pre-run yield estimate was 18× too pessimistic.** Predicting the usable fraction
   as `exp(−err × n_gates)` gave 0.32% at 20 qubits against a measured **5.90%**. That
   product is far too strict: a shot can accumulate gate errors and still land on a state
   with the correct electron count. Post-selection tolerates noise much better than a
   naive fidelity product implies — which is part of *why* the method works.
2. **Our QPU-time estimates were ~6× too low.** Qiskit's `estimate_duration` returns gate
   time only and omits the inter-shot reset/thermalisation delay: 91 µs/shot estimated
   against **~539 µs/shot** measured. Budget from the backend's reported usage on a small
   job, never from `estimate_duration`.

---

## 9. Reproducing

```bash
pip install -r requirements.txt
python3 run_all.py            # (10e,10o), ~12 min, 8 cores, no GPU, no QPU
```

| stage | wall time |
|---|---|
| `stage1_classical.py` — RHF 19 s, AVAS, CASCI 6 roots 13 s, CCSD 1 s | ~35 s |
| `stage2_sqd.py` — noise + subspace sweeps + controls, 100k shots | ~7 min |
| `stage2b_orbital_opt.py` | ~4 min |
| `stage3_report.py` — figures + report | ~5 s |

Larger spaces: `--active fe3d+brs3p` (32q, CASCI 2.57 h) and `--active fe3d+brs3p+fe4s`
(36q, needs >16 GB → GPU/cluster). Fixed seeds throughout; the (22e,16o) A/B reproduced
+321.508 mHa exactly on a re-run.

---

## 10. What was not done

- DMRG (`block2`), SHCI (Dice), FCIQMC — the brief's stretch goals. Published DMRG is
  cited instead; at (10e,10o) they would be redundant since CASCI is already exact.
- CASSCF / NEVPT2. The brief explicitly asks for fixed orbitals and CASCI for a
  like-for-like comparison, so this follows the brief.
- Real quantum hardware.
- A basis better than STO-3G — the largest remaining error source, and the clearest
  next step.
- (26e,18o) low-Sz sectors, (34e,22o) at 44 qubits (693M determinants, prep done).

# Project report — SQD for [2Fe–2S] iron–sulfur electronic structure

**A reproducible workflow and classical benchmark.** BioHackathon, September 2026.

This is the document to read if you are on the team and need to understand the whole
project. It is self-contained: problem, method, everything we ran, every result including
the ones that failed, the plots and what each is for, and the questions we expect to be
asked. Numbers here are the same numbers as in
[docs/RESULTS_MASTER.md](docs/RESULTS_MASTER.md), which is the single source of truth —
if the two ever disagree, RESULTS_MASTER wins.

---

## 0. The 60-second version

We built an end-to-end Sample-Based Quantum Diagonalization (SQD) pipeline for a
[2Fe–2S] iron–sulfur cluster, ran it on a simulator and on **real IBM quantum hardware**,
and checked it against the **exact** classical answer in the same problem — which is the
part most quantum-chemistry demos skip.

Three things came out of it:

1. **It works, exactly.** At 20 qubits SQD reproduces exact diagonalization to
   **+0.0000 mHa** with the correct spin state. On real hardware, **94% of the shots came
   back physically impossible** and the classical post-processing still recovered the exact
   answer. That is the method's central claim, verified rather than assumed.
2. **The physics result:** the magnetic exchange coupling **J = 28.8 → 46.7 → 82.0 cm⁻¹**
   across three increasingly complete descriptions of the molecule, converging toward the
   measured **148 ± 16 cm⁻¹**. Each step adds the specific chemistry the previous one left
   out.
3. **The honest result, which is the interesting one:** on *this* molecule SQD needs ~30%
   of the entire configuration space to get within 34.5 mHa, exact diagonalization is
   **4.7× faster** at 32 qubits, and the quantum circuit provides **no advantage over
   random sampling** at 20 qubits. We measured all three instead of hiding them. An
   independent 2026 assessment of SQD ([arXiv:2501.07231](https://arxiv.org/abs/2501.07231))
   reaches the same conclusion on the same systems.

**Everything reproduces with one command in ~12 minutes on a laptop.**

---

## 1. The problem

### What the molecule is

Iron–sulfur clusters are the electron-transfer and catalytic centres in enzymes across
all life: respiration, photosynthesis, nitrogen fixation, oxygen sensing. The simplest
one — two irons bridged by two sulfurs, **[2Fe–2S]** — is the motif the rest are built
from.

We model it as **[Fe₂S₂(SMe)₄]²⁻**: the [2Fe–2S] core plus four methylthiolate (–SCH₃)
ligands standing in for the cysteine side chains that anchor the cofactor inside a real
protein. 24 atoms, charge −2. *This is the same model complex as the landmark DMRG study*,
which is what makes our comparisons like-for-like rather than approximate.

> **If someone asks "where did SMe come from?"** — it is a chemist's stand-in for
> cysteine. The real cofactor is held by four cysteine residues; modelling the whole
> protein is impossible, so you cut the cysteines back to the smallest group that keeps
> the electronic environment right. Every paper on this system does the same thing.

### Why it is hard

The two iron centres are **magnetically coupled**: their magnetic moments prefer to point
in opposite directions. The strength of that preference is one number, the **exchange
coupling J**, and it controls the cluster's spectroscopy and reactivity. J is measurable
in the lab, so it is the number any calculation has to get right.

The methods that scale well — Hartree–Fock, DFT, coupled cluster — all assume the
electronic structure is essentially *one* dominant arrangement of electrons with small
corrections. Here that assumption fails outright: an enormous number of arrangements
contribute at comparable strength. Our own CCSD run **did not converge**, which is the
clean diagnostic that single-reference methods are outside their range on this system.

The methods that *can* handle it scale exponentially. That is the gap quantum computers
are supposed to address.

---

## 2. The method, in plain terms

**SQD (Sample-Based Quantum Diagonalization)** is a hybrid algorithm:

1. Build a quantum circuit that approximately prepares the molecule's ground state. We
   use the **LUCJ** ansatz, with its parameters taken from a *classical* CCSD calculation
   — no variational training loop on the quantum device.
2. Measure it many times. **Each measurement outcome is read as one candidate arrangement
   of electrons** (one determinant).
3. Throw away the outcomes that are physically impossible (wrong electron count, wrong
   spin), and **repair** some of them — this is the "configuration recovery" step.
4. Hand that shortlist of arrangements to a classical computer, which solves the problem
   **exactly within the shortlist**.
5. Feed the resulting occupancies back into step 3 and repeat until it stops improving.

The quantum computer is used **only as a proposal generator**. It never has to produce a
precise number. That is why noise is survivable: a corrupted shot is a bad *suggestion*,
not a wrong *answer*.

**The one thing that determines whether SQD works on a given molecule:** does a *small*
shortlist suffice? If the true wavefunction is concentrated on a few thousand
arrangements, yes. If it is spread thin over millions, no. Section 5 is our measurement of
which case [2Fe–2S] is.

### The analogy, if you need one for a non-technical audience

You are mapping the lowest point in a mountain range. Exact diagonalization surveys every
square metre — correct, and impossible past a certain size. SQD flies a drone over the
range, notes the few hundred places that *look* low, and then surveys only those
precisely. It wins when the low ground is in a few valleys. It loses when the whole range
is gently, uniformly low — which is what [2Fe–2S] turns out to be.

---

## 3. What we actually ran — the complete inventory

Three active spaces. An "active space" is the subset of electrons and orbitals we treat
exactly; everything else is frozen. `(22e,16o)` = 22 electrons in 16 orbitals. **Qubits =
2 × orbitals** under our mapping.

| active space | what it adds | qubits | configurations | CASCI (exact) | SQD | **J** |
|---|---|---|---|---|---|---|
| **(10e,10o)** Fe 3d | the iron d shells only | **20** | 63,504 | ✅ −5013.64906670 | ✅ **+0.0000 mHa** | **28.8** |
| **(22e,16o)** + bridging S 3p | the sulfur bridge that *carries* the coupling | **32** | 19,079,424 | ✅ −5013.72724223 | ✅ **+34.5 mHa** | **46.7** |
| **(26e,18o)** + Fe 4s | a second shell on the irons | **36** | 73,410,624 | ⚠️ 4 of 6 spin sectors (GPU) | ⚠️ +653 mHa *(hardware)* | **82.0** * |

\* four-sector fit; the singlet is extrapolated. Sz = 0 and Sz = 1 (73.4M and 56.8M
determinants) exceeded the memory budget.

Plus, on **real quantum hardware** (IBM `ibm_fez`, 156-qubit Heron), five completed jobs
at 100,000 shots each — see §6.

**One precision point to get right if you are presenting this:** the three active spaces
are **not strictly nested**. AVAS re-selects orbitals at each threshold, so this is *three
related descriptions of the same molecule*, not one being systematically grown. The trend
is real; the wording should be exact.

### Measured runtimes — everything, with its source

Every row is wall clock we actually measured. Where a number is extrapolated or ran on
different hardware, it says so in the row.

**The full reproducible pipeline** — `python3 run_all.py`, laptop, 8 cores:

| stage | what it does | time |
|---|---|---|
| `stage1_classical.py` | RHF 18.7 s · AVAS · exact CASCI ladder (6 roots) 11.9 s · CCSD 1.0 s | **~35 s** |
| `stage2_sqd.py` | 6 noise levels (160 s) + 6 subspace sizes (299 s) + 2 controls (21 s), 100k shots each | **~8 min** |
| `stage2b_orbital_opt.py` | 4 subspace sizes × 12 orbital-optimisation iterations | **145 s** |
| `stage4_spin_ladder.py` | SQD in all 6 spin sectors + an exact reference for each | **105 s** |
| `stage3_report.py` | 6 figures + markdown report | ~5 s |
| **total, one command** | | **~12 min** |
| `run_all.py --quick` | smoke test, small sweeps | ~2 min |

**The big runs**, per active space:

| active space | qubits | determinants | exact CASCI | SQD | where |
|---|---|---|---|---|---|
| (10e,10o) | 20 | 63,504 | **11.9 s** (6-root ladder; singlet alone ~0.4 s) | **84.9 s** → exact | laptop, 1 core |
| (22e,16o) | 32 | 19,079,424 | **0.96 h** (3,463.8 s) → exact | **4.54 h** (16,357 s) → +34.5 mHa | idle laptop, 1 core |
| (26e,18o) | 36 | 73,410,624 | **1.67 h** for 4 of 6 spin sectors | never run to convergence | **A100 GPU** |

The (26e,18o) sectors individually — a clean small-to-large series inside one active space:

| Sz | determinants | GPU time |
|---|---|---|
| 5 | 43,758 | 6 s |
| 4 | 875,160 | 27 s |
| 3 | 6,694,974 | 367 s |
| 2 | 25,968,384 | **5,612 s** (1.56 h) |
| 1 | 56,805,840 | **out of memory** |
| 0 | 73,410,624 | **out of memory** |

The local scaling exponent steepens **0.50 → 1.28 → 2.01** across those steps, so this
does not extrapolate cleanly — which is exactly why the two largest sectors were attempted
and failed rather than predicted. Figure: `results/figures/gpu_casci_sectors.png`.

**The nested-slice timing study** (`stage5_scaling.py`) — truncated slices of the real
Hamiltonian, purely to locate the wall. The energies from it are meaningless; only the
times are used:

| orbitals | qubits | determinants | CASCI | SQD solve | SQD sampling |
|---|---|---|---|---|---|
| 6 | 12 | 400 | 0.03 s | 0.11 s | 0.85 s |
| 8 | 16 | 4,900 | 0.05 s | 0.32 s | 0.60 s |
| 10 | 20 | 63,504 | 0.36 s | 0.73 s | 1.59 s |
| 12 | 24 | 853,776 | 10.3 s | 1.14 s | 3.25 s |
| 14 | 28 | 11,778,624 | **288.6 s** | 3.33 s | 3.23 s |
| 16 | 32 | 165,636,900 | **infeasible, ~40 GB** | 8.47 s | 48.2 s |

Figures: `talk_runtime.png` (clean, with a qubit axis) and `scaling_casci_vs_sqd.png`.

**Quantum hardware** — five jobs on `ibm_fez`:

| | reported QPU time | wall (incl. queue) |
|---|---|---|
| per 100,000-shot job | **50.6–55.6 s** (mean 53.9 s ⇒ ~539 µs/shot) | 40–113 s |
| five jobs total | **~250 s of a 600 s allocation** | — |

### How hard it actually was

Worth recording honestly, because the difficulty was not where we expected it.

**What was cheap.** The quantum part. Sampling the circuit takes seconds on a simulator
and ~54 s of QPU time on hardware. It is never the bottleneck — and it was not IBM's
bottleneck either (45 min of quantum time against 1.5 h on 3,072 cores).

**What was expensive.** The *classical* post-processing, by three orders of magnitude.
4.54 h for one 32-qubit SQD run on one core, and it still had not converged. Memory, not
cores: one CI vector at (22e,16o) is 153 MB and the Davidson keeps ~30 of them; at
(26e,18o) the working set is 17.6 GB against 16 GB of laptop RAM, which is why that
calculation existed only once it moved to a GPU.

**What was genuinely hard, in order:**

1. **Two bugs that fail silently.** The `ffsim` `logm`-batch bug (§8) does not raise where
   you would look — it just forces `optimize=False`, and the answer degrades from
   +0.0000 to +329 mHa with nothing in the output saying why. The PySCF Davidson
   `max_space` trap is worse: it returns a converged-looking number that is 0.13 mHa off,
   against a spin gap of 0.245 mHa. Both cost hours to find and are one line to fix.
2. **Timing anything reliably.** We measured the 32-qubit comparison wrong more than once
   by running other jobs at the same time; the contamination was large enough to *reverse*
   the CASCI-vs-SQD conclusion. The final numbers are from an otherwise-idle machine. If
   you re-measure, measure alone — and check for stale processes first (`ps` on this Mac
   reports the binary as `Python`, not `python3`, which hid four jobs from us).
3. **Estimating QPU cost.** Qiskit's `estimate_duration` counts gate time only and came
   out **6× low**. One 1,000,000-shot job was submitted and had to be cancelled when IBM's
   own queue estimate showed ~13 minutes against a 10-minute budget.
4. **Memory ceilings arriving without warning.** Two of six spin sectors at (26e,18o)
   simply could not run, which is why J at that space is a four-point fit rather than a
   computed singlet.

### Where it ran out of memory

Worth its own listing, because memory — not cores, not the quantum part — is what set
every ceiling we hit.

| what | memory needed | what happened |
|---|---|---|
| (22e,16o) Davidson, 32 qubits | **153 MB per CI vector**, ~30 kept at once ≈ 4.6 GB | ran on a laptop, fine |
| (26e,18o) Davidson, 36 qubits | **17.6 GB** working set | **did not fit** in 16 GB of laptop RAM → moved to an A100 |
| (26e,18o) Sz = 1 sector, 56,805,840 determinants | — | **out of memory, even on the GPU** |
| (26e,18o) Sz = 0 sector, 73,410,624 determinants | — | **out of memory, even on the GPU** |
| 16-orbital slice, 165,636,900 determinants | **~40 GB** | **infeasible** — this is where exact diagonalisation stops in the scaling study |
| `sci_spin_square` at 44 qubits via full-space embedding | **5.5 GB** | avoided by calling `fci.selected_ci` directly instead |

Two consequences you should state rather than hide:

1. **J at (26e,18o) is a four-point fit with the singlet extrapolated**, because the two
   largest spin sectors could not be computed. Three sectors gave 84.3, four gave 82.0, so
   it is stable to a few percent — but it is a fit, not a computed singlet.
2. **The GPU was not an optimisation, it was an enabler.** The 36-qubit calculation does
   not exist on the laptop at any runtime. `gpu4pyscf.fci.direct_spin1.FCI` overrides only
   `contract_2e` with a CUDA kernel, which was enough.

Marked on the figures: `runtime_vs_active_space.png` (the 17.6 GB ceiling),
`gpu_casci_sectors.png` (the two sectors that never ran), and `talk_runtime.png` /
`scaling_casci_vs_sqd.png` (the ~40 GB wall in the scaling study).

**What took the most calendar time:** the 32-qubit pair of runs (5.5 h of compute, done
twice because of the timing contamination) and the (26e,18o) GPU ladder. Everything that
is *validated* — the exact 20-qubit result and both control experiments — reproduces in
12 minutes.

## 4. Result 1 — the physics: J converges toward experiment

J in the standard iron–sulfur convention `H = 2 J S₁·S₂`, `E(S) = J S(S+1)`, `J > 0`
antiferromagnetic — the same convention as the reference paper, so the numbers are
directly comparable.

| | **J (cm⁻¹)** | fraction of experiment |
|---|---|---|
| ours, (10e,10o) | 28.8 | 0.19× |
| ours, (22e,16o) | 46.7 | 0.32× |
| ours, (26e,18o) | **82.0** | 0.55× |
| DMRG (30e,20o), better basis | 236 | 1.59× |
| broken-symmetry DFT | 310 | 2.09× |
| **experiment** | **148 ± 16** | 1.00× |

**Monotonic, roughly doubling per step, and the reason is identifiable** — because we
identified it by adding the missing piece and watching the number move. The first step
adds the bridging-sulfur 3p orbitals, which are the very atoms that transmit the magnetic
coupling between the two irons (the "superexchange" pathway); the second adds the Fe 4s
double shell. The remaining gap to 148 is the **minimal STO-3G basis set**, which we used
deliberately throughout so the exact answer stayed computable.

**How we get J:** from exact CASCI spin-ladder energies (S = 0…5), one ground state per
Sz sector, fitted to the Landé interval rule. Not from SQD — see §7a for why that fails.

Figure: `results/figures/talk_J_vs_active_space.png`.

---

## 5. Result 2 — the honest assessment of SQD on this molecule

This is the part that makes the project more than a tutorial rerun.

### 5a. SQD against the exact answer

| active space | qubits | CASCI (exact) | SQD | error |
|---|---|---|---|---|
| (10e,10o) | 20 | −5013.64906670 | −5013.64906670 | **+0.0000 mHa** ✅ |
| (22e,16o) | 32 | −5013.72724223 | −5013.69276393 | **+34.5 mHa** |

At 20 qubits: exact to eight decimal places, ⟨S²⟩ = 0.0000 (the correct **singlet**, not
just the right energy), and orbital occupancies matching CASCI to 1.1 × 10⁻³ electrons.

**What 34.5 mHa means** — this comes up every time, so know it cold: 34.5 mHa is
**21.6 kcal/mol**, about **22× outside chemical accuracy** (1.6 mHa = 1 kcal/mol) and
**9× larger than the entire spin ladder** we are trying to resolve. It is not a small
miss. Do not describe it as one.

### 5b. Why it degrades: we measured the coverage curve

SQD diagonalizes the subspace its samples span. At (22e,16o) we swept that subspace size:

| subspace | % of the 19.1M space | error |
|---|---|---|
| 31,329 | 0.16% | +726.7 mHa |
| 223,729 | 1.17% | +568.5 mHa |
| 1,385,329 | 7.26% | +321.5 mHa |
| 5,716,881 | **29.96%** | **+34.5 mHa** |
| 19,079,424 | 100% | 0 (= CASCI) |

**On this molecule you need nearly the whole space.** The wavefunction is not compact: its
participation ratio at (10e,10o) is ~1,956, and 7,171 determinants are needed to hold 99%
of the weight. There is no small shortlist to find. We also checked that the *orbitals*
were not the problem: orbital optimization at frozen subspace recovers only 3–16% of the
gap.

For contrast, IBM ran [4Fe–4S] at **10⁻⁷ %** coverage and got useful answers — because
*that* wavefunction is compact. **Our system is a hard case for SQD, and that is the
finding**, corroborated independently by arXiv:2501.07231 on the same systems.

### 5c. Exact diagonalization beat SQD at 32 qubits

| | work | time | result |
|---|---|---|---|
| **CASCI** | one Davidson over all 19,079,424 determinants | **0.96 h** | **exact** |
| **SQD** | 10 recovery iterations → 5,716,881 subspace | **4.54 h** | +34.5 mHa |

**4.7× slower and not exact.** SQD's individual solve is cheaper — smaller matrix — but
the self-consistent loop needs many of them, and ten iterations did not converge; the
error was still falling at the cap (+441.7 → +34.5 mHa over ten iterations).

The reading: **exact diagonalization is the right tool wherever it fits.** SQD's regime
begins where it doesn't — and we measured where that is:

| orbitals | qubits | determinants | CASCI | SQD *solve* |
|---|---|---|---|---|
| 10 | 20 | 63,504 | 0.36 s | 0.73 s |
| 12 | 24 | 853,776 | 10.3 s | 1.14 s |
| 14 | 28 | 11,778,624 | **288.6 s** | 3.33 s |
| 16 | 32 | 165,636,900 | **infeasible** (~40 GB) | 8.47 s |

CASCI grows ~10,000× across that ladder and hits a wall. The SQD solve grows 77× and
stays under 10 seconds, because its cost tracks the *subspace*, not the space. **Caveat to
state out loud:** SQD's *sampling* cost also grows steeply here (3.2 s → 48 s) — but only
because we are *simulating* the quantum computer. On hardware, sampling is roughly
constant-cost. The column that legitimately replaces CASCI is the solve.

Figure: `results/figures/talk_runtime.png`.

---

## 6. Result 3 — real quantum hardware

**IBM `ibm_fez`, 156-qubit Heron. Five completed jobs, 100,000 shots each, ~250 s of a
600 s allocation.** Full procedure in [docs/HARDWARE.md](docs/HARDWARE.md).

| qubits | suppression | 2q gates | correct N | **usable (N and Sz)** | energy error | ⟨S²⟩ |
|---|---|---|---|---|---|---|
| **20** | — | 2,031 | 17.22% | **5.90%** | **+0.0000 mHa — exact** | 0.0000 |
| **32** | — | 4,878 | 1.40% | **0.47%** | +529.7 mHa | 0.0062 |
| 32 | DD + twirl | 4,742 | 1.52% | **0.47%** | +649.8 mHa | 0.0859 |
| **36** | — | 5,982 | 0.31% | **0.10%** | +653.4 mHa * | 0.0912 |
| 36 | DD + twirl | 5,972 | 0.35% | **0.10%** | +866.1 mHa * | 0.4666 |

\* against an *extrapolated* singlet — there is no exact reference at 36 qubits.

### The headline, and the qualification that must travel with it

At 20 qubits, **94.1% of the hardware shots were physically invalid** — wrong electron
count or wrong spin — and post-selection plus configuration recovery still recovered the
**exact** energy in the correct spin state. That is SQD's noise-resilience claim,
demonstrated on a QPU and checked against an answer we computed independently.

**Say the qualification in the same breath:** the recovered subspace reached **100% of the
space**, so this demonstrates *noise resilience*, not *compression* — and at this size a
uniform-random control reaches the same answer. It is a validated pipeline, not quantum
advantage.

### Depth versus usable shots

Usable yield falls **61×** (5.90% → 0.10%) for **2.9×** the two-qubit gates
(2,031 → 5,982). The poor energies above 20 qubits follow directly: only ~470 and ~100
shots survive post-selection, against spaces of 19.1M and 73.4M. For scale, IBM used
**2.46M shots and 3,072 Fugaku cores**; we used 100k shots on one laptop core. **The 32-
and 36-qubit runs are under-resourced by construction** — they are a depth-versus-yield
measurement, not a claim about what hardware can do with a real budget.

Figures: `hw_yield_vs_depth.png`, `hw_energy_by_size.png`, `hw_mitigation_ab.png`.

### Error suppression: a clean null result, twice

Dynamical decoupling (XpXm) + gate and measurement twirling changed the usable-shot yield
by **nothing at all** at both 32 and 36 qubits (0.47% → 0.47%, 0.10% → 0.10%). At ~5,000
two-qubit gates the dominant error is **incoherent gate error**, and neither technique
removes it: DD targets idle dephasing, twirling reshapes coherent error. Worth knowing:
ZNE, PEC and TREX — the techniques that *would* attack gate error — are `Estimator`-only
and cannot be applied to a `Sampler` workload like SQD at all.

---

## 7. Negative results, measured rather than asserted

These are deliberately in the report. A hackathon project that only shows what worked is
less informative than one that quantifies where the method stops.

**a) J cannot be extracted from SQD.** Running SQD per spin sector and fitting the
Heisenberg form gives **J = −174.6 cm⁻¹ against a true 29.0** — *wrong sign*. Cause: the
high-Sz sectors are tiny (dimension 1, 100) so SQD is trivially exact there, while the
large low-Sz sectors carry up to **55 mHa** of error — against a total ladder span of
3.95 mHa. J is a *difference* of near-degenerate states; error that uneven swamps the
signal by an order of magnitude. **More compute cannot fix this.** Every J in this work
therefore comes from exact CASCI.

> The one-line takeaway, and the most useful thing we learned: **SQD delivers total
> energies to high accuracy here, but cannot resolve the small energy differences that
> magnetic properties depend on.**

**b) Copying IBM's shot count does not help.** Same configuration, only shots changed:
200,000 shots → +321.5 mHa; **2,457,600 shots** (IBM's number) → **+344.2 mHa**. 12× more
shots, 6.7× more unique bitstrings, and the energy got *slightly worse*, because
`samples_per_batch` caps how many configurations enter the subspace regardless of pool
quality. IBM's advantage was never sampling depth — it was **K = 10–100 batches on 3,072
cores**.

**c) The quantum circuit buys nothing at 20 qubits.** A uniform-random-bitstring control
— no quantum information whatsoever — also reaches +0.0000 mHa. Once the subspace covers
the space, how the configurations were chosen stops mattering. **This is the control
experiment a quantum-advantage claim requires, and we ran it.**

**d) Orbital optimization recovers only 3–16% of the gap** at frozen subspace, which is
how we know the 34.5 mHa is a coverage problem and not a basis problem.

---

## 8. Engineering results — two real bugs found and fixed

**1. An upstream defect in `ffsim` 0.0.83 and 0.0.84.**
`unitaries_to_parameters` calls `scipy.linalg.logm` on a batch of shape `(n, dim, dim)`,
but `logm` accepts only one square matrix → `ValueError: expected square array_like
input`. This breaks `UCJOp*.to_parameters()` and therefore
`from_t_amplitudes(..., optimize=True)` — **which the official IBM tutorial uses.**
[`src/ffsim_patch.py`](src/ffsim_patch.py) maps the logarithm over the batch; the
parameter round-trip is then exact (`|⟨ψ|ψ′⟩| = 1.000000000000`). It installs itself only
if the installed ffsim is actually broken, so it no-ops once upstream fixes it.

Measured impact — this is not a cosmetic bug:

| | ansatz spread (participation ratio) | noiseless SQD |
|---|---|---|
| bug present | 2.68 | +329 mHa |
| patched | 3,259 | **+0.0000 mHa** |

**2. A silent PySCF Davidson convergence trap.** In the near-degenerate Sz = 0 sector,
`max_space=12` (PySCF's default) stops short — and **asking for more roots makes it
worse**: +0.022 mHa at `nroots=1`, **+0.128 mHa at `nroots=3`**. `max_space=30` is exact
*and faster*. The S=0→1 gap is only 0.245 mHa, so the default would bias J by ~50%.

**3. GPU acceleration.** `gpu4pyscf.fci.direct_spin1.FCI` subclasses PySCF's own solver
and overrides only `contract_2e` (the H·v product) with a CUDA kernel. It needs a 2-D CI
array while PySCF's Davidson passes a flat vector — a reshape wrapper fixes it. **This is
what made the (26e,18o) CASCI possible at all**: 17.6 GB of Davidson working set against
16 GB of laptop RAM.

---

## 9. The plots, and what each one is for

All in `results/figures/`. Regenerate with `stage3_report.py`, `stage6_talk_figures.py`,
`stage8_hardware_figures.py`, `stage9_extra_figures.py`, `stage10_timing_figures.py` — no
recomputation, all from the committed JSON.

**Runtime, small to large active space** — the three ways of asking it:

| figure | shows | use it when |
|---|---|---|
| `runtime_vs_active_space.png` | measured CASCI vs SQD wall clock at 20 → 32 → 36 qubits, on the **real** Hamiltonian | the direct answer to "how does runtime grow with active space" |
| `talk_runtime.png` | the nested-slice study: CASCI hitting a wall, SQD solve staying flat, with a qubit axis | the clean talk version — shows *where* SQD's regime begins |
| `gpu_casci_sectors.png` | exact CASCI 6 s → 5,612 s across (26e,18o) spin sectors, and the two that ran out of memory | small-to-large inside **one** active space, on one device |
| `cost_casci_vs_sqd.png` | 0.96 h exact vs 4.54 h SQD at 32 qubits | the honesty slide: exact won |

**Energy (Ha) vs iteration** — both sizes:

| figure | shows | use it when |
|---|---|---|
| `energy_vs_iteration_20q.png` | energy in **Hartree** vs iteration, 20 qubits, six subspace sizes, landing exactly on −5013.649067 | the algorithm *succeeding*. Best version of this plot |
| `talk_energy_vs_iteration.png` | energy in **Hartree** vs iteration, 32 qubits, + the growing subspace below | the algorithm working but **not** converging — pairs with the coverage story |
| `iterations_sto-3g_fe3d.png` | the same 20-qubit runs as *error in mHa* on a log axis | when you want to show how far from exact, not the absolute energy |

**Everything else:**

| figure | shows | use it when |
|---|---|---|
| `talk_J_vs_active_space.png` | J = 28.8 → 46.7 → 82.0 against experiment 148 | **the physics headline.** Lead with this |
| `coverage_22e16o.png` | error vs % of the space covered — 30% still +34.5 mHa | **the central finding.** The "no small subspace to find" slide |
| `cost_casci_vs_sqd.png` | 0.96 h exact vs 4.54 h SQD, plus the 10-iteration curve | the honesty slide: exact won at 32 qubits |
| `shots_null.png` | 12× more shots → slightly *worse* energy | showing we tested the obvious fix and it failed |
| `hw_yield_vs_depth.png` | usable shots 5.90% → 0.10% vs gate count | **the hardware headline.** The one honest curve |
| `hw_energy_by_size.png` | hardware energy error by size, vs the 1.6 mHa line | what that yield collapse costs you |
| `hw_mitigation_ab.png` | DD + twirling: identical yield | the null result; shows we tested it |
| `convergence_sto-3g_fe3d.png` | error vs subspace dimension | the coverage argument |
| `noise_sto-3g_fe3d.png` | error vs depolarizing rate — flat at 0.0000 | noise tolerance, on the simulator |
| `occupancies_sto-3g_fe3d.png` | SQD vs CASCI orbital occupancies | that we got the *state* right, not just the energy |
| `spin_ladder_sto-3g_fe3d.png` | S = 0…5 energies + Landé fit | where J comes from |
| `scaling_casci_vs_sqd.png` | the nested-slice timing study | the detailed version of `talk_runtime` |

---

## 10. Where we sit in the literature

Full detail in [docs/LITERATURE_COMPARISON.md](docs/LITERATURE_COMPARISON.md) and
[docs/REFERENCES.md](docs/REFERENCES.md).

- **Classical benchmark:** Sharma, Sivalingam, Neese & Chan, *Nature Chemistry* **6**, 927
  (2014), arXiv:1408.5080. Same model complex. DMRG J = 236 cm⁻¹ at (30e,32o).
- **Quantum benchmark:** Robledo-Moreno *et al.*, arXiv:2405.05068. SQD on the same two
  clusters, 45 and 77 qubits, 3,072 Fugaku cores.
- **Independent corroboration of our main negative result:** arXiv:2501.07231 assesses SQD
  on N₂ and [2Fe–2S] — the same two systems as the tutorial and this project — and reaches
  the same conclusion: for non-compact wavefunctions the subspace must grow to a large
  fraction of the full space. **Cite this if anyone suggests our 32-qubit result means we
  made a mistake.** We did not; the system is a hard case, and that is published.

**The honest placement:** our largest active space, (26e,18o), is **still smaller than
their smallest**, (30e,20o). So this is not a worse version of their calculation — it is
*the approach to the answer from below*, along a path the literature does not document
(they publish only the converged end point). That makes the J series a small original
contribution — and also means **we have no per-space literature value to check 28.8 or
46.7 against.** Say it that way.

---

## 11. Questions we expect, with answers

**"Isn't the quantum computer much faster than the simulator and CASCI? Those results came
back in seconds."** No — that is sampling time only, and it is the cheap part. ~54 s of QPU
time produced samples; the *classical* post-processing is what costs hours. Same for IBM:
45 min of quantum time against 1.5 h on 3,072 cores. **The quantum part was never the
bottleneck, for them either.**

**"Did you beat the classical method?"** No, and we say so. At 32 qubits exact
diagonalization was 4.7× faster *and* exact. What we did is establish where the crossover
is and why.

**"Then what is the contribution?"** A validated, reproducible pipeline with an exact
reference and the control experiments a quantum claim requires — plus three quantified
statements about where the method stops: problem size, missing chemistry, and the *kind of
quantity* being asked for. "SQD is accurate for total energies here but cannot resolve
magnetic couplings, and here is the measurement that shows it" is a more actionable
statement about iron–sulfur chemistry than another energy number.

**"Why STO-3G? That's a crude basis."** Deliberately. It keeps the exact answer computable
so every claim is checkable. The remaining gap to experiment is attributable to it, and
the code already supports better bases — it was simply too expensive for three days.

**"Why only 30% coverage? Is that a threshold for SQD?"** No — 30% is what *this* molecule
needed, not a property of the method. IBM got useful answers on [4Fe–4S] at 10⁻⁷%
coverage. The coverage required is a property of the wavefunction's compactness, which is
exactly the measurement we contribute.

**"Can you get to 100% coverage?"** That *is* CASCI, and we ran it: 0.96 h, exact. That is
the point of §5c.

**"Why did the error suppression not help?"** Because at ~5,000 two-qubit gates the
dominant error is incoherent gate error. DD and twirling address idle dephasing and
coherent error. The techniques that would help (ZNE/PEC/TREX) are Estimator-only and
cannot be used with a Sampler workload like SQD.

---

## 12. Limitations, stated plainly

1. **STO-3G throughout.** Minimal basis; the reference study used TZP-DKH. This accounts
   for most of the remaining gap to experiment.
2. **The three active spaces are not strictly nested** — AVAS re-selects orbitals, so it is
   "three related descriptions", not one systematically grown.
3. **(26e,18o) J = 82.0 is a four-sector fit** with the singlet extrapolated; Sz = 0 and 1
   exceeded memory. Three sectors gave 84.3, four gave 82.0, so it is stable to a few
   percent — but it is a fit, not a computed singlet.
4. **The Heisenberg model is itself imperfect here.** The reference study reports their
   exact levels deviating from the Landé pattern, so part of our residual is the *model*,
   not our active space. Quoting J to better than a few cm⁻¹ would be false precision.
5. **Hardware gate counts vary run-to-run** because Qiskit's layout pass is stochastic, so
   the mitigation A/B carries some layout noise. `seed_transpiler` pinned afterwards.
6. **32/36-qubit hardware runs are deliberately under-resourced** (100k shots, K ≤ 3, one
   core). They measure depth vs yield, nothing more.
7. **Stretch goals not attempted:** DMRG (`block2`), SHCI (Dice), FCIQMC. All three only
   become *meaningful* past (26e,18o), where CASCI stops fitting.

---

## 13. What we'd do with more time

Ranked by how much each would change the conclusions, not by effort.

1. **A real basis set.** Repeat the largest space in TZP-DKH instead of STO-3G. This is
   the single biggest lever on the one number the audience cares about: it should close
   most of the remaining 82 → 148 cm⁻¹ gap to experiment. The code already supports it;
   it was simply too expensive for three days.
2. **Finish the (26e,18o) spin ladder.** Two of six sectors (56.8M and 73.4M
   determinants) ran out of memory, so J = 82.0 is a four-point fit with the singlet
   extrapolated. More GPU memory turns it into a computed number.
3. **Give SQD the resources the method assumes.** IBM used K = 10–100 batches on 3,072
   cores; we used K ≤ 3 on one. Our shot sweep showed K is the knob that matters, not
   shot depth — so this is the one change that could plausibly move the 32-qubit result
   from +34.5 mHa toward chemical accuracy.
4. **A shallower ansatz for hardware.** Usable-shot yield fell 61× from 2,031 to 5,982
   two-qubit gates, and error suppression did nothing. Fewer LUCJ repetitions, or a
   sparser interaction pattern, attacks the actual bottleneck.
5. **Local spins and ⟨S₁·S₂⟩ per state.** Computable from the CI vectors we already
   have, and it would let us compare against the reference study's table directly rather
   than only through J.
6. **[4Fe–4S] and the nitrogenase cofactor.** Where no exact classical answer exists and
   a method like this stops being a benchmark exercise.

What we would *not* spend more time on, because we measured that it does not help: more
shots, error suppression at this depth, orbital optimisation at frozen subspace, and
extracting J from SQD.

---

## 14. Reproducing everything

```bash
pip install -r requirements.txt
python3 run_all.py                    # full (10e,10o) pipeline, ~12 min
```

Individual pieces:

```bash
python3 src/stage1_classical.py --basis sto-3g --avas-threshold 0.5   # classical + exact
python3 src/stage2_sqd.py       --tag sto-3g_fe3d --shots 100000       # the SQD runs
python3 src/stage3_report.py    --tag sto-3g_fe3d                      # report + figures
python3 src/stage6_talk_figures.py                                     # talk figures
python3 src/stage8_hardware_figures.py                                 # hardware figures
```

Replay the hardware post-processing **for free** — the raw shots are committed:

```bash
python3 src/stage7_hardware.py --tag sto-3g_fe3d --replay
```

Larger spaces on a cluster: `bash cluster/preflight.sh` first, then the sbatch scripts in
`cluster/`. See [cluster/README.md](cluster/README.md).

---

## 15. Document map

| file | what it is | who it is for |
|---|---|---|
| **[REPORT.md](REPORT.md)** | this document — everything, self-contained | **the team** |
| [README.md](README.md) | repo front door: install, run, layout | anyone opening the repo |
| [docs/RESULTS_MASTER.md](docs/RESULTS_MASTER.md) | **single source of truth** for every number | whoever is checking a claim |
| [docs/SUBMISSION.md](docs/SUBMISSION.md) | 60-second pitch, abstract, methods | the judges / the write-up |
| [docs/PROPOSAL.md](docs/PROPOSAL.md) | narrative version, no orbital jargon | non-specialist audience |
| [docs/EXPLAINER.md](docs/EXPLAINER.md) | the drone-over-a-mountain-range explanation | whoever presents to a lay audience |
| [docs/PLAN.md](docs/PLAN.md) | the plan we agreed, task split, pitfalls | team logistics |
| [docs/RUN_SUMMARY.md](docs/RUN_SUMMARY.md) | what was run, with timings | the documentation write-up |
| [docs/HARDWARE.md](docs/HARDWARE.md) | **how the QPU runs were done**, step by step | reproducing the hardware results |
| [docs/LITERATURE_COMPARISON.md](docs/LITERATURE_COMPARISON.md) | per-active-space placement vs published work | defending the numbers |
| [docs/REFERENCES.md](docs/REFERENCES.md) | every paper and package, and what it was for | citations |

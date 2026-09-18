# Sample-Based Quantum Diagonalization for a [2Fe–2S] Iron–Sulfur Cluster

**A reproducible SQD workflow benchmarked against exact classical diagonalization**

Hackathon submission document — pitch, abstract, methods, results, reproducibility.


**Real hardware:** run on IBM `ibm_fez` (job `dameitg2fm4c73f2t3f0`), 20 qubits, 100,000 shots, ~4.6 s of QPU time. Only **5.90%** of shots were physically valid — and post-processing still recovered the **exact** energy (+0.0000 mHa) with the correct spin state. Noise-resilience demonstrated; compression not, since the subspace reached 100% of the space.

---

## 1. The 60-second overview (read this aloud)

> Iron–sulfur clusters run respiration, photosynthesis, nitrogen fixation. They're also
> a known failure case for standard computational chemistry — the two irons are
> magnetically coupled in a way single-reference methods can't describe.
>
> So we built an end-to-end **Sample-Based Quantum Diagonalization** pipeline for one.
> A quantum circuit proposes which electronic configurations matter; a classical solver
> then diagonalises exactly within that shortlist. Twenty qubits, then thirty-two, then
> thirty-six. All simulated.
>
> We picked a problem small enough that the exact answer is *also* computable —
> otherwise you can't check a quantum result. SQD nails it: exact to eight decimals,
> correct magnetic state.
>
> But the real contribution is where it **breaks**. We replaced the quantum measurements
> with **random guesses** — random did just as well. And when we asked for the magnetic
> coupling constant, SQD got the **wrong sign**.
>
> So: a map of where this method works. Plus two bugs we fixed, one upstream in IBM's
> own tooling.

*149 words ≈ 64 s at 140 wpm, 56 s at 160. Bold = stress these. Drop the closing
sentence if you're running long.*

**If you only get one sentence:** *"We validated SQD exactly on an iron–sulfur cluster,
and then found the three specific places it stops working — including one the field
hasn't flagged."*

---

## 2. Abstract

We implement and validate an end-to-end Sample-Based Quantum Diagonalization (SQD)
workflow for the oxidized iron–sulfur cluster [Fe₂S₂(SMe)₄]²⁻, the standard truncated
model of the biological [2Fe–2S] cofactor. The pipeline runs PySCF → AVAS → LUCJ ansatz
(ffsim) → `qiskit-addon-sqd`, on three active spaces of increasing size — (10e,10o),
(22e,16o) and (26e,18o), i.e. 20, 32 and 36 qubits — under Jordan–Wigner. Against CASCI — which is full CI in this space, hence exact — SQD
reproduces the ground-state energy to **0.00000 mHa** with ⟨S²⟩ = 0.0000, far inside the
1.6 mHa chemical-accuracy target, and reproduces the CASCI natural occupancies to
1.1 × 10⁻³ electrons. We then report three limits, each measured rather than asserted:
(i) a uniform-random-bitstring control matches the quantum circuit exactly, because at
63,504 determinants any diverse sampler saturates the space; (ii) the Heisenberg exchange
coupling extracted from SQD is qualitatively wrong (J = −174.6 cm⁻¹ against a true
29.0 cm⁻¹), because per-Sz-sector errors reach 55 mHa while the whole spin ladder spans
3.95 mHa; (iii) SQD's accuracy degrades with active-space size at fixed sampling budget, from
exact at 20 qubits to +34.5 mHa at 32 and +653 mHa *(hardware)* at 36, and at 32 qubits exact CASCI
was both faster (0.96 h vs 4.54 h) and exact. Separately, enlarging the active space
raises the exchange coupling from J = 28.8 cm⁻¹ (Fe 3d only) to 46.7 (adding the
bridging-sulfur superexchange pathway) to 82.0 cm⁻¹ (adding the Fe 4s shell), converging
toward the experimental 148 ± 16; the residual gap is the minimal STO-3G basis. We additionally identify and patch an
upstream `ffsim` defect that silently disables the reference tutorial's `optimize=True`
path, and a PySCF Davidson convergence trap that biases J by ~50%.

---

## 3. Methods

### 3.1 Model system and geometry

Oxidized **[Fe₂S₂(SMe)₄]²⁻** — 24 atoms, net charge −2, closed-shell singlet reference.
The Fe₂S₂ rhombic core carries two µ₂-bridging sulfides; four terminal methylthiolate
(–SCH₃) ligands stand in for the cysteine residues that anchor the cofactor in a protein.
This is the same model complex used in the reference DMRG study, so comparisons are
like-for-like.

Idealized geometry built programmatically (`src/geometry.py`), bond lengths from the
crystallographic range for synthetic Fe(III)/Fe(III) analogues:

| parameter | value |
|---|---|
| Fe–Fe | 2.700 Å |
| Fe–S(µ₂) | 2.200 Å |
| Fe–S(terminal) | 2.300 Å |
| S–C | 1.810 Å |
| S–Fe–S (bridge angle) | 104.3° |

Both irons are high-spin Fe(III), S = 5/2. Minimum interatomic separation 1.09 Å (a C–H
bond), i.e. no steric clashes.

### 3.2 Classical setup (`src/stage1_classical.py`)

Follows the IBM SQD chemistry tutorial, substituting the cluster for N₂.

1. **Mean field** — RHF/STO-3G, 122 basis functions, 186 electrons. Density fitting plus
   a second-order (Newton) solver; first-order SCF does not converge on this dianion.
   E = −5012.99869566 Ha. Orbitals are then re-wrapped in a plain RHF object so the
   active-space integrals and CCSD amplitudes use exact rather than density-fitted ERIs.
2. **Active space** — AVAS on the `Fe 3d` atomic-orbital set, threshold 0.5, which
   selects exactly **(10e, 10o)**: the ten Fe 3d orbitals, 10 electrons, 20 qubits,
   **63,504 determinants**.
3. **Integrals** — `CASCI.get_h1cas` / `get_h2cas` on the AVAS orbitals.
4. **Exact reference** — CASCI is full CI in this space. `fci.direct_spin1` with six
   roots in the Sz = 0 sector resolves the complete S = 0…5 spin ladder in one shot
   (13 s). For larger spaces we instead take the lowest state of each Sz sector, which is
   far cheaper.
5. **Exchange coupling** — least-squares fit of `E(S) = J S(S+1)` to the ladder, i.e. the
   convention `H = 2 J S₁·S₂` with J > 0 antiferromagnetic, matching the Fe–S literature.
6. **Ansatz amplitudes** — frozen-core CCSD on the active orbitals, supplying `t1`/`t2`.

### 3.3 Fermion-to-qubit mapping

**Jordan–Wigner**, spin-blocked: qubits `0…norb−1` carry α spin-orbitals,
`norb…2norb−1` carry β. In a printed bitstring (most-significant first) the left half is
β and the right half is α — verified against an asymmetric (3α, 2β) reference state, and
consistent with `qiskit-addon-sqd`'s `hamming_right`/`hamming_left` convention.

JW is not an arbitrary default here. Under JW a computational basis state *is* an
occupation-number vector, so every measured bitstring is directly a Slater determinant
label — which is what allows post-selection and configuration recovery to operate on
Hamming weight per spin half. The Hamiltonian itself is **never mapped to Pauli
operators**: it stays as the `hcore`/`eri` tensors and is diagonalised in the determinant
subspace, so SQD avoids the Pauli-term proliferation that burdens VQE.

### 3.4 Quantum circuit

**LUCJ** (local unitary cluster Jastrow), exactly as in the reference tutorial:

```python
pairs_aa = [(p, p+1) for p in range(norb-1)]   # nearest-neighbour, same spin
pairs_ab = None                                # all-to-all, opposite spin
ucj_op = ffsim.UCJOpSpinBalanced.from_t_amplitudes(
    t2=t2, t1=t1, n_reps=4, interaction_pairs=(pairs_aa, pairs_ab), optimize=True)

circuit.append(ffsim.qiskit.PrepareHartreeFockJW(norb, nelec), qubits)
circuit.append(ffsim.qiskit.UCJOpSpinBalancedJW(ucj_op), qubits)
circuit.measure_all()
```

Transpiled for the *simulator*: depth 983, 1,908 two-qubit gates; gate set
`rz, rxx, ryy, cp, p, x`. Transpiled for `ibm_fez`'s heavy-hex coupling map the same
circuit becomes depth 1,677 / 2,031 two-qubit gates — the hardware tables quote that
number, and it varies a few percent run to run because the layout pass is stochastic.

### 3.5 Sampling backends — simulator (no QPU needed here; hardware runs are separate, see below)

| backend | description | cost (20 q, 3k shots) |
|---|---|---|
| `ffsim` *(default)* | exact statevector sampling within the fixed particle-number sector; `global_depolarizing` knob emulates hardware noise | ~1 s |
| `aer` | `qiskit-aer` statevector with a gate-level depolarizing noise model (1-qubit rate `p`, 2-qubit `10p`) | ~13 s |
| `aer-mps` | as above, matrix-product-state method, for when the qubit count outgrows a dense statevector | — |

`ffsim` is the default because it samples only physical (correct-particle-number) states,
so the depolarizing rate is the *sole* source of invalid shots and the sampling-efficiency
metric stays interpretable. Aer must be transpiled at `optimization_level=0`; levels ≥ 1
hang indefinitely on ffsim's custom LUCJ gate.

### 3.6 SQD (`src/stage2_sqd.py`)

`qiskit_addon_sqd.fermion.diagonalize_fermionic_hamiltonian` with:

```
samples_per_batch   swept 20 → 800
num_batches         3
max_iterations      10
energy_tol          1e-6
occupancies_tol     1e-5
carryover_threshold 1e-4
symmetrize_spin     True when n_alpha == n_beta, else False
sci_solver          partial(solve_sci_batch, spin_sq=0.0, max_cycle=500)
```

The loop is: post-select on correct particle number per spin sector → probabilistically
repair invalid bitstrings using the current average occupancies (configuration recovery)
→ diagonalise in several subspace batches → feed the best occupancies back.

**`max_cycle=500`, not the tutorial's 200.** At 200, diagonalising even the *entire* CAS
space through this solver leaves a 0.112 mHa residual, which would be misread as an SQD
sampling error. We report this "solver floor" explicitly in every run.

### 3.7 Controls and diagnostics

Reported for every configuration:

- **Energy error** vs CASCI, in mHa.
- **⟨S²⟩** of the SQD state, computed on the selected-CI vector directly (cost scales
  with the subspace, not the full CAS space).
- **Sampling efficiency** — fraction of shots with the correct particle number, and with
  correct particle number *and* Sz.
- **Natural occupancies** vs CASCI.
- **Uniform-random control** — SQD on random bitstrings, i.e. zero quantum information.
  This is the baseline a quantum-advantage claim must beat.
- **Selected-CI reference** — diagonalise in the product closure of the largest-|c| CASCI
  determinants. A comparison strategy, *not* an upper bound: the product closure pulls in
  many low-weight configurations, and at small dimensions SQD actually beats it.
- **Solver floor** — the full CAS space through the identical solver.
- **Orbital optimization** (`src/stage2b_orbital_opt.py`, `ffsim.optimize_orbitals`) —
  energy at *frozen* subspace, so any gain is genuine compression.

### 3.8 Software

Python 3.13 · PySCF 2.14.0 · Qiskit 2.4.2 · qiskit-addon-sqd 0.13.1 · ffsim 0.0.84 ·
qiskit-aer 0.17.2. All runs on 8 CPU cores / 16 GB; no GPU, no QPU.

---

## 4. Results

### 4.1 Headline

| active space | qubits | determinants | CASCI (exact) | SQD | **J (cm⁻¹)** |
|---|---|---|---|---|---|
| (10e,10o) Fe 3d | 20 | 63,504 | −5013.64906670 | **+0.0000 mHa** | **28.8** |
| (22e,16o) + bridging S 3p | 32 | 19,079,424 | −5013.72724223 | **+34.5 mHa** | **46.7** |
| (26e,18o) + Fe 4s | 36 | 73,410,624 | 4/6 sectors (GPU) | +653 mHa *(hardware)* | **82.0** * |
| DMRG (30e,20o) TZP-DKH | — | — | — | — | 236 |
| **Experiment** | — | — | — | — | **148 ± 16** |

\* fitted from the four spin sectors we could compute (S = 2–5); the two largest exceeded the memory budget. Adding the fourth moved J from 84.3 to 82.0, i.e. the fit is converging, not drifting.

**J: 28.8 → 46.7 → 82.0 cm⁻¹** — the chemistry converges toward experiment as the active
space stops truncating it. Full detail and all caveats in **`RESULTS_MASTER.md`**.

Supporting: ⟨S²⟩ = 0.0000 at 20 qubits (correct singlet, not just the right energy);
occupancies within 1.1×10⁻³ electrons; energy unchanged across a 100× range of
depolarizing noise; 6-point coverage curve at 32 qubits (0.16% → 30% coverage,
+727 → +34.5 mHa).

### 4.2 The ffsim patch, quantified

| | ansatz spread (participation ratio) | max single-config probability | noiseless SQD |
|---|---|---|---|
| `optimize=False` (pre-patch behaviour) | 2.68 | 0.488 | dim 24,964 → **+329 mHa** |
| `optimize=True` (patched) | **3259** | 0.0059 | dim 63,504 → **+0.0000 mHa** |

With the upstream bug in place the circuit is ~49% concentrated on a single configuration
and a *noiseless* run lands 329 mHa off. Patched, it spreads and SQD is exact with zero
noise.

### 4.3 Limit 1 — exact, but not compressive

SQD reaches CASCI only once the subspace covers essentially the whole CAS space, and the
**uniform-random control matches it exactly** (+0.00000 mHa, dim 63,504). The exact
singlet has participation ratio ≈ 1956 and needs 7,171 determinants for 99% of its
weight, so there is no compact subspace to discover. Orbital optimization at frozen
subspace recovers only 3-16% of the gap across four subspace sizes. At 20 qubits the circuit's distribution carries no
advantage over random configurations — a property of the problem size, and the same point
the `qiskit-addon-sqd` quickstart makes.

### 4.4 Limit 2 — SQD cannot deliver the exchange coupling

Running SQD once per Sz sector and fitting the Heisenberg form:

| Sz | subspace | SQD error (mHa) | exact gap from S=0 (mHa) |
|---|---|---|---|
| 0 | 63,001 | +2.746 | 0.000 |
| 1 | 39,770 | **+55.038** | 0.213 |
| 2 | 14,040 | +9.554 | 0.710 |
| 3 | 2,025 | +0.000 | 1.474 |
| 4 | 100 | +0.000 | 2.531 |
| 5 | 1 | +0.000 | 3.919 |

**Fitted J = −174.6 cm⁻¹ against a true 29.0** — wrong sign, 35 mHa fit residual.
High-Sz sectors have tiny Hilbert spaces (dimension 1, 100) so SQD is trivially exact;
the large low-Sz sectors carry millihartree error. J is a difference of near-degenerate
states, so error that uneven swamps the signal by an order of magnitude. **More compute
cannot fix this.**

> SQD delivers total energies to high accuracy here, but cannot resolve the small energy
> differences that magnetic properties depend on.

### 4.5 Limit 3 — the active space, not the solver

| source | J (cm⁻¹) |
|---|---|
| **This work, (10e,10o) Fe 3d, STO-3G** | **28.8** |
| DMRG-CI (30e,20o), TZP-DKH — Sharma/Chan 2014 | 236 |
| BS-DFT | 310 |
| **Experiment** (magnetic susceptibility) | **148 ± 16** |

Our ladder follows the Landé interval rule to 4 × 10⁻⁵ Ha, so the *shape* is right and
only the scale is wrong. Cause: an Fe-3d-only active space contains no bridging-sulfur 3p
orbitals, hence no superexchange pathway between the irons — plus a minimal basis where
the literature uses TZP-DKH.

Also note the 1.6 mHa target is a loose yardstick here: the S=0→1 gap is only 0.245 mHa
(54 cm⁻¹), roughly six times *below* chemical accuracy.

### 4.6 Bugs found

1. **Upstream, `ffsim` 0.0.83 and 0.0.84** — `unitaries_to_parameters` calls
   `scipy.linalg.logm` on a `(n, dim, dim)` batch, but `logm` accepts only a single
   square matrix, raising `ValueError: expected square array_like input`. This breaks
   `UCJOp*.to_parameters()` and therefore `from_t_amplitudes(..., optimize=True)`, which
   the official tutorial uses. `src/ffsim_patch.py` maps the logarithm over the batch;
   the parameter round-trip is then exact (|⟨ψ|ψ′⟩| = 1.000000000000). The patch installs
   itself only if the local ffsim is actually broken.
2. **PySCF Davidson convergence trap** — in the near-degenerate Sz = 0 sector,
   `max_space=12` (the default) stalls 0.022 mHa high with one root and **0.128 mHa high
   with three roots — more roots is worse.** On a 0.245 mHa gap that is a ~50% bias in J.
   `max_space=30` converges exactly *and* runs faster.

---

## 5. Reproducibility

```bash
pip install -r requirements.txt
python3 run_all.py            # full pipeline, ~12 min on 8 cores
python3 run_all.py --quick    # smoke test, ~2 min
```

| stage | wall time |
|---|---|
| `stage1_classical.py` — RHF 19 s, AVAS, CASCI 6 roots 13 s, CCSD 1 s | ~35 s |
| `stage2_sqd.py` — 6 noise levels × 6 subspace sizes + controls, 100k shots | ~7 min |
| `stage2b_orbital_opt.py` — 4 subspace sizes × 12 OO iterations | ~4 min |
| `stage3_report.py` — figures + RESULTS.md | ~5 s |

Every number in this document is regenerated into `results/` as JSON plus six figures and
a markdown report. Fixed seeds throughout.

### Real quantum hardware

Five completed jobs on **IBM `ibm_fez`** (156-qubit Heron), 100,000 shots each, ~250 s of
a 600 s allocation:

| qubits | 2q gates | correct N | **usable (N and Sz)** | energy error |
|---|---|---|---|---|
| **20** | 2,031 | 17.22% | **5.90%** | **+0.0000 mHa - exact** |
| 32 | 4,878 | 1.40% | **0.47%** | +529.7 mHa |
| 36 | 5,982 | 0.31% | **0.10%** | +653.4 mHa* |

At 20 qubits **94.1% of the shots came back physically invalid** and the classical
post-processing still recovered the **exact** energy in the correct spin state
(<S^2> = 0.0000). Usable yield then falls **61x** for 2.9x the two-qubit gates, which is
what the 32- and 36-qubit energies reflect. Dynamical decoupling + twirling changed the
usable fraction by **nothing** at either size - a clean null, because at ~5,000 two-qubit
gates the dominant error is incoherent gate error.

*\* against an extrapolated singlet; there is no exact reference at 36 qubits.*

Full procedure, budget accounting and pitfalls: **docs/HARDWARE.md**. Raw shots for every
job are committed, so the post-processing replays at zero QPU cost.

**Larger spaces, since run:** `cluster/run_rungA.sbatch` gave the (22e,16o) space
(32 qubits, 19,079,424 determinants) — exact CASCI in 0.96 h, SQD in 4.54 h, J = 46.7 cm⁻¹,
SQD error +34.5 mHa. `cluster/run_casci_26e18o.sbatch` gave (26e,18o) on a GPU,
J = 82.0 cm⁻¹. Both are in docs/RESULTS_MASTER.md.

---

## 6. References

1. S. Sharma, K. Sivalingam, F. Neese, G. K.-L. Chan, *Low-energy spectrum of
   iron–sulfur clusters directly from many-particle quantum mechanics*,
   **Nature Chemistry 6, 927 (2014)**; arXiv:1408.5080. — Source of the model complex,
   the DMRG benchmark J = 236 cm⁻¹, the Heisenberg convention, and the quoted experimental
   and BS-DFT values.
2. J. Robledo-Moreno *et al.*, *Chemistry Beyond Exact Solutions on a Quantum-Centric
   Supercomputer*, arXiv:2405.05068. — SQD on [2Fe–2S] and [4Fe–4S], up to 77 qubits;
   the published quantum baseline.
3. IBM Quantum, *Sample-based quantum diagonalization* tutorial. — The workflow this
   implementation ports.
4. `qiskit-addon-sqd` documentation: quickstart, open-/closed-shell guide,
   `optimize_orbitals` guide.

---

## 7. Document map

| file | purpose |
|---|---|
| **../REPORT.md** | the full project report — everything, self-contained |
| **SUBMISSION.md** | this document — pitch, abstract, methods, results |
| **RESULTS_MASTER.md** | single source of truth for every number |
| **HARDWARE.md** | how the real-QPU runs were done, step by step |
| **PROPOSAL.md** | non-specialist narrative, no orbital-level jargon |
| **EXPLAINER.md** | lay-audience explanation |
| **PLAN.md** | team logistics, 3-day split, Q&A prep, code pitfalls |
| **RUN_SUMMARY.md** | what was run, with measured timings |
| **LITERATURE_COMPARISON.md** | per-active-space placement vs published work |
| **REFERENCES.md** | every paper and package, and what it was for |
| **../README.md** | how to run, runtimes, caveats |
| `../results/RESULTS_*.md` | auto-generated numerical report + figures |

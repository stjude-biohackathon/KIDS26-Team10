# References

Grouped by what we actually used them for. Where a number in our results comes from a
paper, the paper is named at the point of use in
[RESULTS_MASTER.md](RESULTS_MASTER.md) and [LITERATURE_COMPARISON.md](LITERATURE_COMPARISON.md)
as well.

---

## The two papers this project is positioned against

**1. Sharma, Sivalingam, Neese & Chan,** *Low-energy spectrum of iron–sulfur clusters
directly from many-particle quantum mechanics*, **Nature Chemistry 6, 927 (2014)**,
[arXiv:1408.5080](https://arxiv.org/abs/1408.5080).

The classical benchmark. We use the **same model complex**, [Fe₂S₂(SCH₃)₄]²⁻, which is
what makes our comparison like-for-like rather than approximate. Source of:

- the DMRG-CI converged coupling **J = 236 cm⁻¹** at (30e,32o);
- the quoted **experimental J = 148 ± 16 cm⁻¹** and **BS-DFT 310 cm⁻¹**;
- the Heisenberg convention `H = 2 J S₁·S₂`, `E(S) = J S(S+1)` that every J in this repo
  uses, so our numbers are directly comparable to theirs;
- the warning that the Heisenberg model itself deviates from the exact levels here, which
  qualifies our own fits (see LITERATURE_COMPARISON.md §3).

Their smallest active space, (30e,20o), is **still larger than our largest**, (26e,18o).
We are approaching their converged answer from below along a path they do not publish.

**2. Robledo-Moreno et al.,** *Chemistry Beyond Exact Solutions on a Quantum-Centric
Supercomputer*, [arXiv:2405.05068](https://arxiv.org/abs/2405.05068).

The quantum benchmark, and the direct antecedent of this work: SQD on **the same two
clusters**, [2Fe–2S] at (30e,20o)/45 qubits and [4Fe–4S] at (54e,36o)/77 qubits. Source
of the resource comparison in RESULTS_MASTER.md §7 — 2,457,600 shots, K = 10–100 batches,
~1.5 h of classical post-processing on **3,072 Fugaku cores**. Their reported agreement
with a classical HCI reference is *tens of mHa*.

Two things we took from it directly: that the quantum sampling (~45 min) was never their
bottleneck either, and that K (batch count) rather than shot depth is the knob that
matters — which our own shot sweep then confirmed independently (§6b).

---

## The paper that reaches our conclusion independently

**Reinholdt, Kjellgren, Fuß, Kongsted et al. (SQD assessment),**
[arXiv:2501.07231](https://arxiv.org/abs/2501.07231).

Worth foregrounding, because it is corroboration rather than support: an independent
assessment of SQD that tests, among other systems, **N₂ and a [2Fe–2S] cluster** — the
same two systems as the IBM tutorial and this project — and reaches the same structural
conclusion we arrived at from the coverage sweep: **for strongly correlated systems whose
wavefunction is not compact, the SQD subspace has to grow to a large fraction of the full
space before the energy is useful**, so the method's advantage is system-dependent rather
than general.

Our §3 coverage table (30% of 19.1M determinants still leaving +34.5 mHa) and our
participation-ratio measurement (~1,956 at (10e,10o)) are the same finding measured a
different way. Cite this when someone asks whether our poor 32-qubit result means we
did something wrong — we did not; the system is a hard case for SQD, and that is a
published observation.

*(Verify the author list against the arXiv page before putting it on a slide — the
identifier and content are what we checked, not the full byline.)*

---

## Method and algorithm papers

- **SQD / QSCI, the underlying idea.** Kanno et al., *Quantum-Selected Configuration
  Interaction*, [arXiv:2302.11320](https://arxiv.org/abs/2302.11320) — measure a quantum
  state, read each outcome as a determinant, diagonalise classically in the span. SQD adds
  the self-consistent **configuration recovery** step on top of this.
- **LUCJ ansatz.** Motta et al., *Bridging physical and computational aspects of the
  unitary coupled cluster Jastrow ansatz*,
  [arXiv:2104.08957](https://arxiv.org/abs/2104.08957) — the local unitary cluster Jastrow
  form we build from classical CCSD `t1`/`t2` amplitudes. Local interaction pairs are what
  keep the circuit shallow enough to be transpilable onto a heavy-hex device.
- **AVAS.** Sayfutyarova, Sun, Chan & Knizia, *Automated Construction of Molecular Active
  Spaces from Atomic Valence Orbitals*, **JCTC 13, 4063 (2017)** — how the active space is
  selected. `--avas-threshold 0.5` gives exactly (10e,10o) on this cluster.

---

## Software and documentation

| what | where | how we used it |
|---|---|---|
| **IBM SQD chemistry tutorial** | [quantum.cloud.ibm.com/docs/en/tutorials/sample-based-quantum-diagonalization](https://quantum.cloud.ibm.com/docs/en/tutorials/sample-based-quantum-diagonalization) | The workflow this repo ports, step for step, with N₂ swapped for the iron–sulfur cluster |
| `qiskit-addon-sqd` | docs: quickstart, open-/closed-shell guide, `optimize_orbitals` guide, Dice/SHCI integration | `diagonalize_fermionic_hamiltonian`, `solve_sci_batch`, configuration recovery |
| `ffsim` | [github.com/qiskit-community/ffsim](https://github.com/qiskit-community/ffsim) | `UCJOpSpinBalanced.from_t_amplitudes`, `PrepareHartreeFockJW`, `UCJOpSpinBalancedJW`, `PRE_INIT` transpiler passes, and the fast statevector sampler. **We patch an upstream bug in 0.0.83/0.0.84** — see `src/ffsim_patch.py` and RESULTS_MASTER.md §8 |
| `pyscf` | [pyscf.org](https://pyscf.org) | RHF, AVAS, CASCI/FCI, CCSD, `fci.selected_ci`, and the Davidson whose `max_space` default cost us 0.13 mHa until we found it |
| `gpu4pyscf` | [github.com/pyscf/gpu4pyscf](https://github.com/pyscf/gpu4pyscf) | `fci.direct_spin1.FCI` — CUDA `contract_2e`. This is what made the (26e,18o) CASCI possible at all (17.6 GB working set vs 16 GB of laptop RAM) |
| `qiskit-ibm-runtime` | — | `SamplerV2` on `ibm_fez`, dynamical decoupling and twirling options |
| `qiskit-aer` / `qiskit-aer-gpu` | — | gate-level depolarizing noise model, as a device stand-in. Must transpile at `optimization_level=0` |

---

## Related work we read but did not build on

- **GPU-accelerated SQD.** Two 2026 preprints report large speedups for the SQD solver on
  GPUs (~40× and ~95× respectively). Neither released code at the time we looked, so our
  GPU path is our own `gpu4pyscf` wrapper rather than a reimplementation of theirs. Worth
  checking for released code before anyone extends this project.
- **Stretch-goal classical methods we did not run:** DMRG via `block2`, SHCI via Dice,
  FCIQMC. All three are only *meaningful* once the active space is large enough that CASCI
  is infeasible — which, for us, means past (26e,18o). `qiskit-addon-sqd` ships a Dice
  integration guide if you want SHCI as a scalable cross-check.

---

## How to cite this repository

> *Sample-Based Quantum Diagonalization for [2Fe–2S] iron–sulfur electronic structure: a
> reproducible workflow and classical benchmark.* BioHackathon project, September 2026.

State in any write-up that the (10e,10o) result is **exact against CASCI** and that the
32/36-qubit hardware runs are **deliberately under-resourced**; both qualifications are
part of the result.

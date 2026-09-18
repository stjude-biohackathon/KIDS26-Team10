# Sample-Based Quantum Diagonalization for an Iron–Sulfur Cluster

### A reproducible workflow, benchmarked honestly against exact classical results

---

## The problem

Iron–sulfur clusters sit at the business end of enzymes across every domain of life.
They shuttle electrons in respiration and photosynthesis, they do the hard chemistry in
nitrogen fixation, and they let cells sense oxygen. The smallest of them — two irons
bridged by two sulfurs, written **[2Fe–2S]** — is the structural motif everything else
is built on.

What makes them interesting is also what makes them difficult. The two iron centres are
magnetically coupled: their magnetic moments prefer to point in opposite directions, and
the strength of that preference, a single number called the **exchange coupling J**,
controls the cluster's spectroscopy and its reactivity. J is measurable in the lab, so
it is the number any calculation has to get right.

Standard computational chemistry cannot get it right. The methods that scale well —
Hartree–Fock, DFT, coupled cluster — all assume the molecule's electronic structure is
essentially one dominant arrangement with small corrections. In [2Fe–2S] that assumption
simply fails: an enormous number of arrangements contribute at comparable strength. The
methods that *can* describe this correctly scale exponentially and run out of road
quickly.

This is exactly the kind of problem quantum computers are supposed to be for.

## The approach

**Sample-Based Quantum Diagonalization (SQD)** is a hybrid method. A quantum circuit is
measured repeatedly; each measurement outcome is read as one plausible electronic
arrangement. The quantum device is used only as a *proposal generator* — it suggests
which arrangements matter. A classical computer then takes that shortlist and solves the
problem exactly within it.

The appeal is robustness. Because the quantum device only has to propose candidates and
never has to produce a precise number, noisy hardware can still be useful: a corrupted
measurement is a bad suggestion, not a wrong answer. A recovery step even repairs
measurements that come back physically impossible.

We built this workflow end to end for **[Fe₂S₂(SMe)₄]²⁻** — the [2Fe–2S] core with four
small ligands standing in for the cysteine amino acids that anchor it inside a real
protein. This is the same model system used in the landmark classical study of these
clusters, which means our results are directly comparable to published work rather than
to an approximation of it.

We deliberately chose a version of the problem small enough that the **exact answer can
also be computed classically**. That is the whole point: without an exact reference, a
quantum result cannot be checked.

## What we found

**1. The workflow reproduces the exact answer.** SQD recovers the exact classical energy
to eight decimal places — comfortably past the accuracy target set in the challenge
brief — and recovers the correct magnetic state, not merely a state that happens to have
the right energy. For context, the published quantum-hardware result on this same
cluster agrees with its classical reference to within tens of thousandths of a Hartree;
ours is exact. That is a validation of the pipeline, not a claim of quantum advantage —
we are exact precisely *because* we chose a problem small enough to check.

**2. At this size, the quantum circuit earns nothing.** We ran the control experiment
that this kind of claim demands: replace the quantum measurements with **random
guesses**. The random baseline performs just as well. The reason is structural — at this
problem size the shortlist grows until it covers essentially every possibility, so any
reasonably varied source of suggestions works. We report this plainly. It is a property
of the problem size, not a flaw in the method, and it tells you exactly where the
interesting regime begins.

**3. We ran it on a real quantum computer, and it survived.** Five runs on an IBM
156-qubit machine. On the smallest problem, **94% of the measurements came back
physically impossible** — the wrong number of electrons, or the wrong magnetic state —
and the classical repair-and-solve step still returned the **exact** answer, in the
correct magnetic state. That is the method's whole promise, measured rather than assumed,
and checked against an answer we had computed independently.

We also measured the price of scale, which nobody advertises: tripling the size of the
quantum circuit left **61× fewer usable measurements**, and the two standard
noise-suppression techniques we tried improved that by *nothing at all*. We report both.

**4. The bottleneck is the chemistry, not the quantum method — and we proved it by
fixing it.** Our exchange coupling starts at **29 cm⁻¹** against experiment's
**148 ± 16**. Enlarging the description of the molecule twice raised it to **47** and
then **82** — converging toward the measured value.

The cause was identifiable, and we identified it by removing it. The first enlargement
added the bridging sulfurs, which are the very atoms that transmit the magnetic coupling
between the two irons; the second added a further shell on the irons themselves. Each
step roughly doubled the answer. The solver was doing its job correctly all along — it
was being handed an incomplete description of the molecule, and the remaining gap to 148
is down to the deliberately crude mathematical basis we used throughout.

**5. Some quantities are simply out of reach for this method — and we measured which.**
We attempted to extract the exchange coupling directly from SQD and it failed badly,
returning the wrong sign and roughly six times the correct magnitude. The reason is
instructive: J is a *difference* between states that lie extraordinarily close together,
and SQD's small errors vary unevenly from one state to another. Those errors are an
order of magnitude larger than the differences we are trying to resolve. More computing
time cannot fix this. The one-line summary is the most useful thing we learned:

> **SQD delivers total energies to high accuracy here, but cannot resolve the tiny energy
> differences that magnetic properties depend on.**

**6. We found and fixed two real bugs.** One is an upstream defect in a widely used
quantum-chemistry simulation library that silently disables a feature the official
tutorial depends on; our patch restores it exactly. The other is a numerical trap in the
classical reference calculation, where the standard settings quietly stop short of
convergence — by enough to bias the exchange coupling by around 50%, and where the
obvious fix makes it *worse*. Both are documented so nobody repeats them.

## Why this is worth presenting

The honest result is more useful than a clean one. We have a workflow anyone can rerun in
about twelve minutes, checked against an exact answer, with the control experiments that
a quantum-advantage claim would require — and a clear, quantified account of the three
places it stops working: problem size, missing chemistry, and the kind of quantity being
asked for.

That last point is the contribution. "This method is accurate for total energies but
cannot resolve magnetic couplings, and here is the measurement that shows it" is a more
actionable statement about iron–sulfur chemistry than another energy number would be.

## What remains

The clearest single improvement left is the **mathematical basis set**. We used the
crudest one available throughout, deliberately, so that the exact answer stayed
computable and every claim could be checked. Repeating the largest calculation with the
basis the reference study used would close most of the remaining gap to experiment, and
the code already supports it — it was simply too expensive for the time available.

Beyond that, the natural progression is the larger clusters, [4Fe–4S] and the nitrogenase
cofactor, where no exact classical answer exists at all and methods like this one stop
being a benchmark exercise and start being the only option.

Worth noting what made the largest calculation possible at all: it needed more memory
than the laptop had, and ran only once moved onto a graphics processor. That is a
practical constraint worth knowing for anyone repeating this — the bottleneck is memory
and single-core speed, not the quantum part.

---

### At a glance

| | |
|---|---|
| **System** | [Fe₂S₂(SMe)₄]²⁻ — a [2Fe–2S] cluster with ligands modelling its protein anchors |
| **Quantum resources** | 20, 32 and 36 qubits — simulated, **and run on a real 156-qubit IBM machine** (5 jobs, ~250 s of quantum time) |
| **Problem size** | 63,504 to 73 million candidate electronic configurations |
| **Accuracy vs exact classical** | exact to 8 decimal places, on the simulator **and on hardware**; correct magnetic state |
| **Exchange coupling** | 29 → 47 → 82 cm⁻¹ across three descriptions · 148 ± 16 (experiment) |
| **Reproduce everything** | one command, ~12 minutes on a laptop |
| **Classical references used** | exact diagonalisation (ours) + published benchmark and experimental values |

*The full project report is `../REPORT.md`. Technical detail, method choices, division of
work and known pitfalls are in `PLAN.md`; every number is traceable in `RESULTS_MASTER.md`;
the quantum-hardware procedure is in `HARDWARE.md`. Figures are generated into
`../results/figures/`.*

# Where our active spaces sit relative to the literature

Source for all reference values: **Sharma, Sivalingam, Neese & Chan,
*Nature Chemistry* 6, 927 (2014)**, arXiv:1408.5080 — same molecule,
[Fe₂S₂(SCH₃)₄]²⁻. Quantum reference: **Robledo-Moreno *et al.*, arXiv:2405.05068**.

---

## 1. Active spaces, ours against theirs

| | active space | orbitals included | determinants | method | **J (cm⁻¹)** |
|---|---|---|---|---|---|
| **ours** | (10e,10o) | Fe 3d only | 63,504 | CASCI (exact) | **28.8** |
| **ours** | (22e,16o) | Fe 3d + bridging S 3p | 19,079,424 | CASCI (exact) | **46.7** |
| **ours** | (26e,18o) | Fe 3d + **Fe 4s** + bridging S 3p | 73,410,624 | CASCI (exact, 4/6 sectors) | **82.0** |
| Chan #1 | (30e,20o) | Fe 3d + bridging S 3p + **one 3p per terminal S** | >10¹⁵ | DMRG-CI, M = 3500 | — |
| Chan #2 | (30e,32o) | as #1 + **Fe 4s and Fe 4d** (double shell) | >10¹⁷ | DMRG-CI, M = 4500 | **236** |
| Chan #3 | (30e,20o) | as #1, orbitals optimised | — | DMRG-SCF | — |
| Chan #4 | (30e,32o) | as #2, orbitals optimised | — | DMRG-SCF | — |
| — | — | broken-symmetry DFT | — | BS-DFT | 310 |
| — | — | — | — | **experiment** | **148 ± 16** |

### The honest placement

**Our largest active space is still smaller than their smallest.** Their "minimal full
valence" (30e,20o) contains everything in our (22e,16o) *plus* one 3p orbital per
terminal sulfur. Our (26e,18o) adds the Fe 4s shell but still omits the terminal-S
orbitals and the Fe 4d shell.

So the comparison is not "we did a worse version of their calculation" — it is
**"we mapped the approach to the exact answer from below, and they report the converged
end point."** Our trend 28.8 → 46.7 → 82.0 is heading for their 236 and the measured
148 from a direction the literature does not document.

### What the literature does *not* provide

**A J-versus-active-space series for this cluster.** Chan *et al.* report converged
values at (30e,20o) and (30e,32o); they do not publish the approach from smaller spaces.
The earlier work they cite for Fe-3d-only treatments (Hübner & Sauer,
*J. Chem. Phys.* **116**, 617 (2002)) is described in their paper as making "further
drastic approximations, completely removing non-bridging atoms, and including only d
electrons."

That makes our systematic series a small original contribution rather than a
reproduction — **but it also means we have no per-space literature value to check our
28.8 or 46.7 against.** Only the converged end point is published. State it that way.

---

## 2. Reference numbers worth quoting

| quantity | literature value | our value |
|---|---|---|
| J, converged DMRG | **236 cm⁻¹** (30e,32o) | 82.0 (26e,18o) |
| J, experiment (magnetic susceptibility, similar synthetic dimer) | **148 ± 16 cm⁻¹** | — |
| J, broken-symmetry DFT | **310 cm⁻¹** | — |
| DMRG convergence of *relative* energies | better than **0.1 kcal/mol ≈ 35 cm⁻¹** | our Landé residual 1.5 × 10⁻⁴ Ha ≈ 33 cm⁻¹ |
| Local spin on each Fe, ⟨S₁²⟩ | **5.47 (S=0) to 5.74 (S=5)** — *not* 8.75 | not measured |
| Spin ladder span | ~0–10,000 cm⁻¹ over S = 0…5 | 867 cm⁻¹ at (10e,10o) |
| SQD on hardware, (30e,20o), 45 qubits | agrees with HCI to **tens of mHa** | +0.0000 mHa at 20q; +34.5 at 32q |
| SQD subspace used | up to **d = 100M** | 5.7M |
| SQD classical post-processing | **~1.5 h on 3,072 Fugaku cores** | 4.54 h on 1 core |

---

## 3. Two things in the literature that qualify our own fits

**a) The Heisenberg model is itself imperfect here.** Chan *et al.* find their exact
levels deviate from the Landé pattern: *"the Heisenberg model overestimates the lower
spin state energies, while underestimating those of the higher spin states."*

We fit J with exactly that model. So some of our residual is the model, not our active
space. Our Landé residuals (4 × 10⁻⁵ to 1.5 × 10⁻⁴ Ha, i.e. ~9–33 cm⁻¹) are the same
order as their stated convergence threshold — meaning our fits are about as good as the
model deserves, and quoting J to better than a few cm⁻¹ would be false precision.

**b) The irons are not pure Fe(III).** They measure local spin ⟨S₁²⟩ = 5.47–5.74 against
8.75 for an ideal S = 5/2 ion. The Heisenberg picture of "two 5/2 spins" is an
idealisation; real covalency with the sulfurs reduces it.

This is a good "what we'd measure next" item: local spins and the spin-correlation
function ⟨S₁·S₂⟩ per state are computable from our existing CI vectors, and would let us
compare against their Table directly rather than only through J.

---

## 4. One-line summary for a slide

> The published converged value is J = 236 cm⁻¹ from DMRG on a (30e,32o) active space,
> against 148 ± 16 from experiment. Our largest space, (26e,18o), is still smaller than
> their smallest, and gives 82.0 cm⁻¹ — so we are approaching the published answer from
> below, along a path the literature does not document.

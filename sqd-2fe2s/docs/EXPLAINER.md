# Explaining this to people who don't do quantum chemistry

Use this for the talk. No jargon, in the order that makes it land.

---

## 1. The molecule (15 seconds)

Iron–sulfur clusters are tiny metal clumps that sit inside enzymes. They run
respiration, photosynthesis, and nitrogen fixation. Life is full of them.

Ours is the smallest: **two iron atoms bridged by two sulfurs.**

The two irons act like two tiny magnets, and they prefer to point *opposite* to each
other. **How strongly they prefer that is one number, called J** — and you can measure
it in a lab. For our cluster, experiment says **148**.

So we have a scoreboard: compute J, compare to 148. That's the whole problem.

---

## 2. Why it's hard (30 seconds)

To get J you need to know how the electrons arrange themselves.

Normal chemistry software assumes there's **one main arrangement** plus small
corrections. That works for most molecules. It completely fails here — in this cluster
a huge number of arrangements matter *equally*.

**We didn't just assert that, we showed it:** the standard workhorse method (coupled
cluster) doesn't even converge on this molecule. It gives up.

The methods that *do* work have to consider every possible arrangement at once — and the
number of arrangements explodes:

| how much of the molecule you describe | arrangements to consider |
|---|---|
| a little | 63,504 |
| more | 19,079,424 |
| more still | 73,410,624 |
| the amount you'd actually want | more than there are stars in the galaxy |

That explosion is why people want a quantum computer here.

---

## 3. What SQD does — use this analogy

> **Finding the lowest point in a mountain range.**
>
> **The exact method** surveys every square metre. Guaranteed right, but the area doubles
> every time you widen the map, so you run out of time and memory.
>
> **SQD** flies a drone over the range first. The drone isn't precise — it just points at
> the valleys worth surveying. Then you survey *only those valleys*, exactly.
>
> The quantum computer is the drone. It doesn't compute the answer; it **suggests where
> to look.** All the precision comes from the exact survey afterwards, on a normal
> computer.

Why that's clever: a drone with a blurry camera is still useful. **A noisy quantum
computer still gives useful suggestions** — a bad suggestion is just a wasted survey, not
a wrong answer. That's why this method works on today's imperfect hardware.

---

## 4. The catch, and this is our actual result

The drone trick only works **if the lowest point really is in a few valleys.**

If instead the range is one enormous shallow basin, the drone has to point at nearly
everything — and you've saved nothing.

**Our cluster turned out to be the shallow basin.** We measured it:

| fraction of arrangements searched | how wrong the answer was |
|---|---|
| 0.16% | 727 |
| 7% | 322 |
| **30%** | **34** |
| 100% | 0 (this *is* the exact method) |

*(units: millihartree — "chemical accuracy", the threshold for a useful answer, is 1.6)*

So even searching **30%** of the arrangements, we were still **22× short** of a useful
answer. To actually get it right, we needed all of it — at which point we were just doing
the exact method the slow way.

**For contrast:** the IBM team ran a bigger cluster searching **0.0000001%** of its
arrangements and got useful answers — because *their* system was the few-deep-valleys
kind.

### The one-line version

> **The method needs the answer to be concentrated in a few places. We found a case
> where it isn't, and measured exactly what that costs.**

That is a genuinely useful finding. "Here is when this works and when it doesn't" is more
actionable than one more energy number.

---

## 5. The part that did work beautifully

Remember the scoreboard: experiment says J = **148**.

| how much of the molecule we described | J we computed |
|---|---|
| just the iron atoms | **29** |
| + the sulfur bridge between them | **47** |
| + one more layer on the irons | **82** |
| **experiment** | **148** |

Each time we described more of the molecule, the answer roughly **doubled toward the
right value.** The sulfur bridge mattered because it is physically the thing that
transmits the magnetism between the two irons — take it out of the model, and the model
can't know the irons talk to each other.

**The maths was never wrong. The model was incomplete, and we proved it by completing it.**

The remaining gap to 148 is one more known approximation we ran out of time to remove.

---

## 6. Questions you will get

**"Did the quantum computer help?"**
Honestly, at the size we could check: no. We ran the control — replaced the quantum
suggestions with **random guesses**, and got the same answer. At this problem size
anything varied enough works. That control is the thing that makes the rest credible, and
most people don't run it.

**"So why bother?"**
Because of the runtime plot. The exact method hits a wall — it runs out of memory at
166 million arrangements. The subspace method keeps running. Everything we could *check*
was below the wall, by design; the interesting regime is above it.

**"Did you use a real quantum computer?"**
**Yes — five runs on an IBM 156-qubit machine**, and most of the work on a simulator
alongside it. Both on purpose. On hardware you find out whether the method survives real
noise; on a simulator you can turn the noise off and see exactly what the quantum part
contributed, which is the only way to tell whether it mattered at all.

**"What happened on the real machine?"**
On the smallest problem, **94 out of every 100 measurements came back impossible** — the
wrong number of electrons, or the wrong magnetic state. The repair-and-solve step still
returned the **exact** answer. That is the entire promise of the method, and it held.

Then we made the circuit three times bigger and watched the useful measurements drop by a
factor of **61**. We also tried the two standard noise-suppression tricks, and they
changed the number that matters by *nothing*. Both of those are results too.

**"Is 34 millihartree good?"**
No. It's 21 kcal/mol, about a chemical bond. We say so.

---

## 7. Slide order that works

1. **The molecule and the scoreboard** — two irons, J, experiment says 148
2. **Why normal methods fail** — coupled cluster doesn't converge; arrangements explode
3. **The drone analogy** — what SQD is
4. **`talk_runtime.png`** — exact hits a wall, subspace method doesn't
5. **`talk_energy_vs_iteration.png`** — watch it converge on the exact answer
6. **`talk_J_vs_active_space.png`** — 29 → 47 → 82 vs 148. *The payoff slide.*
7. **`hw_yield_vs_depth.png`** — the real machine: 94% of shots junk, exact answer anyway,
   then the 61× collapse with depth. *The credibility slide.*
8. **The catch** — the coverage table; the method needs a concentrated answer
9. **What we'd do next** — remove the last approximation

Lead with the molecule and the scoreboard. Never open with "active space".

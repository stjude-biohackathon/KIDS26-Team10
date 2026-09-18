# How the quantum-hardware runs were done

Everything needed to reproduce the five real-QPU jobs in
[RESULTS_MASTER.md](RESULTS_MASTER.md#real-quantum-hardware--ibm-ibm_fez-156-qubit-heron),
including the mistakes that cost us budget. The code is one file:
[`src/stage7_hardware.py`](../src/stage7_hardware.py).

**Backend:** `ibm_fez`, 156-qubit IBM Heron.
**Budget:** a 10-minute (600 s) allocation on the IBM Quantum Platform Open/paid plan.
**Spent:** ~250 s across five completed jobs.

---

## 1. Credentials — do this yourself, in your own shell

Get the API key and the CRN (instance identifier) from
<https://quantum.cloud.ibm.com>. Then, in a Python shell **on your own machine**:

```python
from qiskit_ibm_runtime import QiskitRuntimeService
QiskitRuntimeService.save_account(
    channel="ibm_quantum_platform",
    token="<YOUR_API_KEY>",
    instance="<YOUR_CRN>",
    name="default",
    set_as_default=True,
    overwrite=True,
)
```

This writes `~/.qiskit/qiskit-ibm.json`. Nothing in this repository contains or reads a
token — `stage7_hardware.py` calls `QiskitRuntimeService()` with no arguments and picks up
the saved account. **Do not paste your key into a chat window, a notebook you will commit,
or a shell one-liner that lands in your history.**

Check it works without spending anything:

```bash
python3 -c "from qiskit_ibm_runtime import QiskitRuntimeService as S; print([b.name for b in S().backends(operational=True)])"
```

---

## 2. What actually gets submitted

Exactly the same circuit the simulator runs, transpiled for the device. From
`stage7_hardware.py`:

```python
pairs_aa = [(p, p + 1) for p in range(norb - 1)]      # nearest-neighbour α–α
ucj = ffsim.UCJOpSpinBalanced.from_t_amplitudes(
    t2=data["t2"], t1=data["t1"], n_reps=n_reps,
    interaction_pairs=(pairs_aa, None), optimize=True)  # needs src/ffsim_patch.py

qr = QuantumRegister(2 * norb, "q")
circuit = QuantumCircuit(qr)
circuit.append(ffsim.qiskit.PrepareHartreeFockJW(norb, nelec), qr)
circuit.append(ffsim.qiskit.UCJOpSpinBalancedJW(ucj), qr)
circuit.measure_all()
```

- **LUCJ ansatz** built from the *classical* CCSD `t1`/`t2` amplitudes computed in
  stage 1 — no variational optimisation on hardware, no parameter training loop. The
  circuit is fixed before submission, which is what makes a single Sampler job enough.
- **Jordan–Wigner, spin-blocked**: qubits `0 … norb-1` are α, `norb … 2·norb-1` are β.
  In a printed Qiskit bitstring the **left half is β and the right half is α**. Getting
  this backwards silently halves your usable-shot count.
- Transpilation uses ffsim's own pass set:

  ```python
  pm = generate_preset_pass_manager(optimization_level=3, backend=backend,
                                    pre_init=ffsim.qiskit.PRE_INIT)
  ```

  `PRE_INIT` is not optional — it merges the orbital-rotation and Givens structure before
  the generic passes see it. Note also that on the **simulator** path Aer must run at
  `optimization_level=0`; levels ≥ 1 hang on ffsim's custom LUCJ gate.

Resulting hardware circuits:

| qubits | active space | transpiled depth | two-qubit gates |
|---|---|---|---|
| 20 | (10e,10o) | 1,677 | 2,031 |
| 32 | (22e,16o) | 3,123 | 4,878 |
| 36 | (26e,18o) | 3,425 | 5,982 |

Gate counts vary a few percent run-to-run because Qiskit's layout pass is stochastic.
Pin `seed_transpiler` if you want them comparable across runs.

---

## 3. The procedure, in the order we ran it

### Step 1 — estimate first. Always.

```bash
python3 src/stage7_hardware.py --tag sto-3g_fe3d --shots 100000 --estimate
```

`--estimate` transpiles, prints depth/gate counts and Qiskit's `estimate_duration`, and
**exits without submitting**. Zero QPU cost.

> **Read this before you trust the number.** `estimate_duration` counts *gate time only*.
> It omits the inter-shot reset and thermalisation delay, which dominates. Ours came out
> **91 µs/shot** against **~539 µs/shot** actually billed — a 6× underestimate. Use
> `--estimate` for gate counts and a sanity check, then calibrate from the `quantum_seconds`
> the backend reports on one small real job.

### Step 2 — submit one job

```bash
python3 src/stage7_hardware.py --tag sto-3g_fe3d --shots 100000
```

One circuit, one `SamplerV2` job, no sweeps. The script prints the job ID immediately,
waits, then prints the reported usage the moment it lands.

The three sizes we ran:

```bash
python3 src/stage7_hardware.py --tag sto-3g_fe3d                  --shots 100000  # 20q
python3 src/stage7_hardware.py --tag "sto-3g_fe3d+brs3p"          --shots 100000  # 32q
python3 src/stage7_hardware.py --tag "sto-3g_fe3d+brs3p+fe4s"     --shots 100000  # 36q
```

Each needs its `results/stage1_<tag>.npz` integrals to exist already — the hardware stage
does no chemistry of its own.

### Step 3 — error suppression A/B

```bash
python3 src/stage7_hardware.py --tag "sto-3g_fe3d+brs3p" --shots 100000 \
        --dd --dd-sequence XpXm --twirling --suffix _mitigated
```

`--suffix` keeps the output files from overwriting the plain arm. Result: a clean null —
see [RESULTS_MASTER.md](RESULTS_MASTER.md). Note what is *not* available: ZNE, PEC and
TREX are `Estimator`-only. SQD is a `Sampler` workload, so they cannot be applied to it at
all, no matter how poor the raw samples look.

### Step 4 — post-process for free, as often as you like

Every job writes its raw shots to `results/hardware_samples_<tag>[<suffix>].json`
**before** any post-processing runs. So all the classical work — post-selection,
configuration recovery, subspace diagonalisation — replays at zero QPU cost:

```bash
python3 src/stage7_hardware.py --tag sto-3g_fe3d --replay
python3 src/stage7_hardware.py --tag "sto-3g_fe3d+brs3p" --replay --suffix _mitigated
```

This is the single most useful design choice in the script. Three of our five results
were produced or corrected by replay after the QPU time was already gone.

### Step 5 — collect and plot

```bash
python3 src/collect_hardware.py          # -> results/hardware_all.json
python3 src/stage8_hardware_figures.py   # -> results/figures/hw_*.png
```

---

## 4. Budget accounting

| job | qubits | suppression | reported `quantum_seconds` | wall |
|---|---|---|---|---|
| `dameitg2fm4c73f2t3f0` | 20 | — | not reported | 40 s |
| `damemt78gn2s739larvg` | 32 | — | not reported | 58 s |
| `dameqg8pqrnc7397es80` | 32 | DD + twirl | **55.64** | 98 s |
| `dameumn8gn2s739lbb70` | 36 | — | **50.64** | 47 s |
| `damf48gpqrnc7397f9t0` | 36 | DD + twirl | **55.45** | 113 s |

Mean **53.9 s per 100,000-shot job** ⇒ **~539 µs/shot**. Estimated total **~250 s of
600 s**. Wall time includes queueing and is not billed.

**One job cancelled:** `damf1k02fm4c73f2trag`, 1,000,000 shots at 32 qubits. IBM's own
queue estimate showed ~13 minutes — over the entire remaining allocation — so it was
cancelled before execution. Watch the platform's estimate after submitting, not just
Qiskit's.

---

## 5. Pitfalls, in the order they will bite you

1. **`estimate_duration` is gate time only.** ~6× low. Calibrate from real reported usage.
2. **Save the raw samples before post-processing.** Post-processing takes minutes to hours
   at 32–36 qubits and can crash on an edge case; losing QPU-derived shots to a local
   `KeyError` is unforgivable.
3. **Check the platform's queue estimate after submit.** Cancel if it exceeds budget.
4. **`--max-full-space-dim` does not exist in stage 7.** It is a stage-2 flag. Passing it
   fails in `argparse` before any submission — harmless, but check your command twice
   anyway, because a typo that *is* a valid flag will spend real money.
5. **Bitstring halves: left = β, right = α.** Verify your post-selection against the known
   exact case (20 qubits) before trusting a size where you have no reference.
6. **Error suppression will not rescue a deep circuit.** DD and twirling address idle
   dephasing and coherent error; at ~5,000 two-qubit gates the dominant term is incoherent
   gate error. Shorten the circuit or use fewer qubits instead.

---

## 6. What the hardware runs do and do not show

**Do:** at 20 qubits, 94.1% of shots came back physically invalid and the classical
post-processing still recovered the **exact** energy (+0.0000 mHa) in the correct spin
state (⟨S²⟩ = 0.0000). That is SQD's noise-resilience claim, measured on a QPU against an
answer we computed independently.

**Do not:** show compression or advantage. The recovered subspace at 20 qubits is 100% of
the space, and a uniform-random-bitstring control reaches the same energy. The 32- and
36-qubit runs are deliberately under-resourced (100k shots, K ≤ 3, one core, against IBM's
2.46M shots and 3,072 cores) and are reported as a depth-versus-yield measurement, not as
a statement about what hardware can achieve with a real budget.

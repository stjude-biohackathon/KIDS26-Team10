"""Stage 7 - run the LUCJ circuit on real quantum hardware, then post-process locally.

Purpose: a genuine QPU data point, on the one active space where we have the **exact**
answer, so the hardware result is checkable rather than decorative.

    (10e,10o) Fe 3d, 20 qubits, exact CASCI singlet = -5013.64906670 Ha

The demo is the central claim of SQD: hardware samples will be badly corrupted (our
circuit is 983 deep with ~1900 two-qubit gates, so most shots come back with the wrong
electron count), and the classical post-processing is supposed to recover a good energy
anyway. Either it does or it doesn't - both are a result.

QPU TIME IS THE SCARCE RESOURCE. This script therefore:
  * submits ONE circuit, never a sweep
  * defaults to 10,000 shots - deliberately small; check the reported usage, then decide
  * prints the usage the moment the job finishes
  * saves the raw samples to disk so the SQD post-processing (which is free and local)
    can be re-run any number of times without touching the QPU again

Usage:
    python3 src/stage7_hardware.py --shots 10000            # submit, ~seconds of QPU
    python3 src/stage7_hardware.py --replay                 # re-post-process, no QPU
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time
from functools import partial

import ffsim
import numpy as np
from qiskit import QuantumCircuit, QuantumRegister

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import ffsim_patch  # noqa: E402,F401
from common import RESULTS, load_json, save_json  # noqa: E402
from stage2_sqd import SCI_MAX_CYCLE, run_sqd, sampling_efficiency  # noqa: E402

DEFAULT_TAG = "sto-3g_fe3d"


def build(data, norb, nelec, n_reps):
    pairs_aa = [(p, p + 1) for p in range(norb - 1)]
    ucj = ffsim.UCJOpSpinBalanced.from_t_amplitudes(
        t2=data["t2"], t1=data["t1"], n_reps=n_reps,
        interaction_pairs=(pairs_aa, None), optimize=True)
    qr = QuantumRegister(2 * norb, "q")
    c = QuantumCircuit(qr)
    c.append(ffsim.qiskit.PrepareHartreeFockJW(norb, nelec), qr)
    c.append(ffsim.qiskit.UCJOpSpinBalancedJW(ucj), qr)
    c.measure_all()
    return c


def transpile_for(backend, circuit):
    """Map the circuit onto the device.  Local only - costs no QPU time."""
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
    pm = generate_preset_pass_manager(optimization_level=2, backend=backend)
    # ffsim ships passes that map LUCJ onto heavy-hex connectivity, which is what the
    # ansatz was designed for.  Two-qubit gates are the fidelity budget, so this matters.
    pm.pre_init = ffsim.qiskit.PRE_INIT
    isa = pm.run(circuit)
    ops = dict(isa.count_ops())
    n2q = sum(v for k, v in ops.items() if k in ("cz", "cx", "ecr", "rzz"))
    return isa, ops, n2q


def estimate(args, circuit, norb, nelec):
    """Print the expected QPU cost WITHOUT submitting anything.

    Transpilation and backend properties are metadata reads; neither consumes the
    time budget.  Run this before spending, always.
    """
    from qiskit_ibm_runtime import QiskitRuntimeService
    service = QiskitRuntimeService()
    backend = (service.backend(args.backend) if args.backend
               else service.least_busy(operational=True, simulator=False,
                                       min_num_qubits=2 * norb))
    isa, ops, n2q = transpile_for(backend, circuit)

    print(f"backend           {backend.name}  ({backend.num_qubits} qubits)")
    print(f"transpiled depth  {isa.depth()}")
    print(f"two-qubit gates   {n2q}")
    print(f"shots requested   {args.shots:,}")

    # per-shot duration from the instruction schedule, if the backend exposes it
    per_shot = None
    try:
        per_shot = isa.estimate_duration(backend.target, unit="s")
    except Exception:
        try:
            dt = backend.configuration().dt
            per_shot = isa.duration * dt if isa.duration else None
        except Exception:
            pass
    if per_shot:
        qpu = per_shot * args.shots
        print(f"\nper-shot circuit  {per_shot * 1e6:.1f} us")
        print(f"est. QPU time     {qpu:.1f} s  ({qpu / 600 * 100:.1f}% of a 600 s budget)")
        print("  (excludes queue wait and per-job overhead, which are not QPU time)")
    else:
        print("\nbackend did not expose a duration estimate.")

    # fidelity expectation, so the result is not a surprise
    try:
        errs = [backend.target[g][q].error for g in ("cz", "ecr", "cx")
                if g in backend.target for q in backend.target[g]
                if backend.target[g][q].error]
        if errs:
            med = sorted(errs)[len(errs) // 2]
            import math
            print(f"\nmedian 2q error   {med:.2e}")
            print(f"expected fraction of shots with no 2q error: "
                  f"{math.exp(-med * n2q) * 100:.2f}%")
            print("  post-selection + configuration recovery is meant to rescue the rest")
    except Exception:
        pass
    print("\nNO QPU TIME USED. Re-run without --estimate to submit.")


def submit(args, circuit, norb, nelec):
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2

    service = QiskitRuntimeService()
    backend = (service.backend(args.backend) if args.backend
               else service.least_busy(operational=True, simulator=False,
                                       min_num_qubits=2 * norb))
    print(f"backend: {backend.name}  ({backend.num_qubits} qubits)", flush=True)
    isa, ops, n2q = transpile_for(backend, circuit)
    print(f"transpiled for hardware: depth {isa.depth()}, {n2q} two-qubit gates",
          flush=True)
    print(f"submitting {args.shots:,} shots ...", flush=True)

    sampler = SamplerV2(mode=backend)

    # Error suppression available to SamplerV2.  Note the headline mitigation methods
    # (ZNE, PEC, TREX) are Estimator-only: they return mitigated *expectation values*,
    # and SQD needs raw bitstrings, so they cannot be used here.
    if args.dd:
        # Dynamical decoupling fills idle windows with pulse sequences that refocus
        # dephasing.  This circuit is 3000+ deep, so many qubits idle a long time.
        sampler.options.dynamical_decoupling.enable = True
        sampler.options.dynamical_decoupling.sequence_type = args.dd_sequence
        print(f"  dynamical decoupling ON ({args.dd_sequence})", flush=True)
    if args.twirling:
        # Twirling randomises coherent error into stochastic Pauli error.  That is worth
        # more here than it sounds: configuration recovery is built to repair random bit
        # flips, so turning coherent drift into random flips plays to its strengths.
        sampler.options.twirling.enable_gates = True
        sampler.options.twirling.enable_measure = True
        print("  gate + measurement twirling ON", flush=True)

    t0 = time.time()
    job = sampler.run([isa], shots=args.shots)
    print(f"job id {job.job_id()} - waiting", flush=True)
    result = job.result()
    wall = time.time() - t0

    # The usage figure is the thing to watch: it is what comes off the 10-minute budget.
    usage = None
    try:
        usage = job.usage_estimation or job.metrics().get("usage", {})
    except Exception:
        pass
    print(f"\nwall {wall:.0f}s   reported usage: {usage}", flush=True)

    bit_array = result[0].data.meas
    save_json(RAW, {
        "backend": backend.name, "job_id": job.job_id(), "shots": args.shots,
        "transpiled_depth": int(isa.depth()), "transpiled_2q": int(n2q),
        "wall_s": wall, "usage": str(usage),
        # store raw bytes so post-processing can be replayed without the QPU
        "array": bit_array.array.tolist(), "num_bits": int(bit_array.num_bits),
    })
    print(f"raw samples saved to {RAW} - post-processing can be replayed for free")
    return bit_array


def load_raw():
    from qiskit.primitives.containers import BitArray
    d = load_json(RAW)
    arr = np.asarray(d["array"], dtype=np.uint8)
    print(f"replaying {d['shots']:,} shots from {d['backend']} (job {d['job_id']})")
    return BitArray(arr, d["num_bits"]), d


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tag", default=DEFAULT_TAG,
                    help="which stage-1 active space to run on hardware")
    ap.add_argument("--shots", type=int, default=10_000,
                    help="keep this small; QPU time is the scarce resource")
    ap.add_argument("--backend", default=None, help="else least busy")
    ap.add_argument("--n-reps", type=int, default=4)
    ap.add_argument("--samples-per-batch", type=int, default=300)
    ap.add_argument("--num-batches", type=int, default=3)
    ap.add_argument("--max-iterations", type=int, default=5)
    ap.add_argument("--dd", action="store_true",
                    help="dynamical decoupling: suppress dephasing on idling qubits")
    ap.add_argument("--dd-sequence", default="XpXm",
                    choices=["XX", "XpXm", "XY4"])
    ap.add_argument("--twirling", action="store_true",
                    help="gate + measurement twirling: turn coherent error into random "
                         "bit flips, which configuration recovery is designed to repair")
    ap.add_argument("--suffix", default="",
                    help="tag the output files, e.g. _mitigated, so a mitigated run does "
                         "not overwrite the unmitigated baseline")
    ap.add_argument("--estimate", action="store_true",
                    help="print the expected QPU cost and expected shot yield, and "
                         "submit nothing. Uses no budget. Do this first.")
    ap.add_argument("--replay", action="store_true",
                    help="re-post-process saved hardware samples; touches no QPU")
    args = ap.parse_args()

    RAW = RESULTS / f"hardware_samples_{args.tag}{args.suffix}.json"
    globals()["RAW"] = RAW
    s1 = load_json(RESULTS / f"stage1_{args.tag}.json")
    data = np.load(RESULTS / f"stage1_{args.tag}.npz")
    hcore, eri = data["hcore"], data["eri"]
    e_nuc = float(data["nuclear_repulsion_energy"])
    norb = s1["norb"]; nelec = tuple(s1["nelec_active"])
    e_exact = s1["e_casci_singlet"]
    ref = (f"exact CASCI = {e_exact:.8f} Ha" if e_exact is not None
           else "NO exact reference (space past exact CI)")
    print(f"{args.tag}: ({sum(nelec)}e,{norb}o)  {2 * norb} qubits | {ref}\n")

    if args.estimate:
        estimate(args, build(data, norb, nelec, args.n_reps), norb, nelec)
        return

    if args.replay:
        bit_array, meta = load_raw()
    else:
        circuit = build(data, norb, nelec, args.n_reps)
        bit_array = submit(args, circuit, norb, nelec)
        meta = load_json(RAW)

    eff = sampling_efficiency(bit_array, norb, nelec)
    print(f"\nRAW HARDWARE SAMPLES")
    print(f"  shots                      {eff['shots']:,}")
    print(f"  correct electron count     {eff['frac_correct_particle_number'] * 100:.2f}%")
    print(f"  correct count AND spin     {eff['frac_correct_particle_number_and_sz'] * 100:.2f}%")
    print(f"  unique bitstrings          {eff['unique_bitstrings']:,}")

    res = run_sqd(hcore, eri, e_nuc, bit_array, norb, nelec,
                  args.samples_per_batch, args.num_batches, args.max_iterations,
                  seed=2026, spin_sq=0.0)
    err = (res["energy"] - e_exact) * 1e3 if e_exact is not None else float("nan")
    print(f"\nAFTER CLASSICAL POST-PROCESSING")
    print(f"  energy                     {res['energy']:.8f} Ha")
    print(f"  error vs exact             {err:+.4f} mHa")
    print(f"  subspace                   {res['subspace_dim']:,} "
          f"({100 * res['subspace_dim'] / s1['cas_dim']:.1f}% of the space)")
    print(f"  <S^2>                      {res['spin_sq']:.4f}")

    save_json(RESULTS / f"stage7_hardware_{args.tag}{args.suffix}.json", {
        "backend": meta["backend"], "job_id": meta["job_id"], "shots": meta["shots"],
        "transpiled_depth": meta["transpiled_depth"],
        "transpiled_2q": meta["transpiled_2q"], "usage": meta.get("usage"),
        "e_casci_exact": e_exact, "raw_efficiency": eff,
        "sqd_energy": res["energy"], "sqd_error_mha": err,
        "sqd_subspace_dim": res["subspace_dim"], "sqd_spin_sq": res["spin_sq"],
        "dynamical_decoupling": args.dd, "twirling": args.twirling,
    })
    print(f"\nwrote results/stage7_hardware_{args.tag}{args.suffix}.json")


if __name__ == "__main__":
    main()

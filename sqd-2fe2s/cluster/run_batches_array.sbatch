#!/bin/bash
#SBATCH --job-name=sqd-batches
#SBATCH --array=0-9
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=24G
#SBATCH --time=08:00:00
#SBATCH --output=batch-%A_%a.out --error=batch-%A_%a.err
#
# =============================================================================
#  K=10 INDEPENDENT BATCHES - copying IBM's setup, the part that needs a cluster
# =============================================================================
#
# IBM (arXiv:2405.05068) used K=10 batches for [2Fe-2S] and K=100 for [4Fe-4S].
# We had been running K=1-3.  Batches are independent subspaces; SQD keeps the
# best.  More batches = better odds of hitting a good subspace.
#
# They are perfectly parallel, and our solver uses exactly ONE core (measured
# 98.7% CPU over a 4 h run).  So ten batches as ten single-core tasks finish in
# the wall time of one - which is precisely what a cluster is for, and what a
# single 32-core job would waste.
#
# Also copies IBM's sampling depth: 2,457,600 shots instead of our 200,000.
#
# Each task uses a different seed, so each explores a different subspace.
# Collect with cluster/collect_batches.py afterwards.
#
# Target: (22e,16o), 32 qubits, which HAS an exact reference
# (-5013.72724223 Ha), so every batch gets a real error bar.
# Baseline to beat: +34.5 mHa (K=3, 200k shots, 10 iterations, 4.54 h on 1 core).
# =============================================================================
set -euo pipefail
module load python/3.13 2>/dev/null || true
source "${VENV:-$HOME/venvs/sqd}/bin/activate" 2>/dev/null || true
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"; mkdir -p results/batches

TAG=sto-3g_fe3d+brs3p
SEED=$((2026 + SLURM_ARRAY_TASK_ID))

python3 -u - <<PY
import sys, time, json; sys.path.insert(0,"src")
import ffsim_patch, numpy as np, ffsim
from functools import partial
from qiskit import QuantumCircuit, QuantumRegister
from qiskit_addon_sqd.fermion import diagonalize_fermionic_hamiltonian, solve_sci_batch

seed = $SEED; task = $SLURM_ARRAY_TASK_ID
E_EXACT = -5013.72724223; CAS = 19079424
d = np.load("results/stage1_${TAG}.npz")
h, eri = d["hcore"], d["eri"]; en = float(d["nuclear_repulsion_energy"])
norb = 16; nelec = (11, 11)
pa = [(p, p+1) for p in range(norb-1)]

ucj = ffsim.UCJOpSpinBalanced.from_t_amplitudes(
    t2=d["t2"], t1=d["t1"], n_reps=4, interaction_pairs=(pa, None), optimize=True)
qr = QuantumRegister(2*norb); c = QuantumCircuit(qr)
c.append(ffsim.qiskit.PrepareHartreeFockJW(norb, nelec), qr)
c.append(ffsim.qiskit.UCJOpSpinBalancedJW(ucj), qr); c.measure_all()

# IBM's sampling depth
ba = ffsim.qiskit.FfsimSampler(default_shots=2_457_600, norb=norb,
        nelec=nelec, seed=seed).run([c]).result()[0].data.meas

t0 = time.time()
r = diagonalize_fermionic_hamiltonian(
    h, eri, ba, samples_per_batch=800, norb=norb, nelec=nelec,
    num_batches=1, max_iterations=5,
    sci_solver=partial(solve_sci_batch, spin_sq=0.0, max_cycle=500),
    symmetrize_spin=True, seed=seed)
e = float(r.energy + en); dim = int(np.prod(r.sci_state.amplitudes.shape))
out = {"task": task, "seed": seed, "shots": 2457600, "samples_per_batch": 800,
       "energy": e, "err_mha": (e - E_EXACT)*1e3, "subspace_dim": dim,
       "pct_cas": 100*dim/CAS, "wall_s": time.time()-t0}
json.dump(out, open(f"results/batches/task{task:02d}.json", "w"), indent=2)
print(f"task {task} seed {seed}: E={e:.8f} err={(e-E_EXACT)*1e3:+.3f} mHa "
      f"dim={dim:,} ({100*dim/CAS:.2f}% CAS) {time.time()-t0:.0f}s")
PY

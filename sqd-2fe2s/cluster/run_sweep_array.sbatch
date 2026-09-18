#!/bin/bash
#SBATCH --job-name=sqd-sweep
#SBATCH --array=0-7
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=sweep-%A_%a.out --error=sweep-%A_%a.err
#
# =============================================================================
#  SQD coverage sweep as a JOB ARRAY - one samples_per_batch value per task
# =============================================================================
#
# Why an array and not one big job: the selected-CI solver is **single-threaded**
# (measured: 98.7% CPU = one core, on a run that lasted 4 h).  Giving one job 32
# cores therefore wastes 31 of them.  Eight independent single-core tasks finish
# a whole sweep in the time one point would take.
#
# This also means a Xeon cluster node is *not* faster than a modern laptop for a
# single point - per-core clock is what matters here.  The cluster wins on
# concurrency and memory, not on speed.
#
#   --cpus-per-task=1   one core is all a task can use
#   --array=0-7         eight coverage points at once
#   --mem=24G           per task; the subspace is the memory driver
#
# Set SPACE below to choose the active space.
# =============================================================================
set -euo pipefail
module load python/3.13 2>/dev/null || true
source "${VENV:-$HOME/venvs/sqd}/bin/activate" 2>/dev/null || true
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1        # keep tasks from fighting
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"; mkdir -p results/sweep

SPACE=${SPACE:-fe3d+brs3p+fe4s}          # (26e,18o) 36 qubits
TAG=sto-3g_${SPACE}
SPBS=(50 100 200 400 800 1600 3200 6400)
SPB=${SPBS[$SLURM_ARRAY_TASK_ID]}

# task 0 also produces the integrals the others need; the rest wait for them
if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ] && [ ! -f "results/stage1_${TAG}.npz" ]; then
    python3 -u src/stage1_classical.py --basis sto-3g --active "$SPACE" \
            --skip-casci --tag "$TAG"
fi
for i in $(seq 1 120); do
    [ -f "results/stage1_${TAG}.npz" ] && break
    echo "task $SLURM_ARRAY_TASK_ID waiting for integrals ($i)"; sleep 15
done

python3 -u - <<PY
import sys, time, json; sys.path.insert(0,"src")
import ffsim_patch, numpy as np, ffsim
from functools import partial
from scipy.special import comb
from qiskit import QuantumCircuit, QuantumRegister
from qiskit_addon_sqd.fermion import diagonalize_fermionic_hamiltonian, solve_sci_batch

spb = int("$SPB")
d = np.load("results/stage1_${TAG}.npz")
h, eri = d["hcore"], d["eri"]; en = float(d["nuclear_repulsion_energy"])
norb = h.shape[0]; ne = int(round(d["occ_a"].sum() + d["occ_b"].sum()))
nelec = (ne//2, ne - ne//2)
CAS = int(comb(norb, nelec[0]))**2
pa = [(p, p+1) for p in range(norb-1)]

ucj = ffsim.UCJOpSpinBalanced.from_t_amplitudes(
    t2=d["t2"], t1=d["t1"], n_reps=4, interaction_pairs=(pa, None), optimize=True)
qr = QuantumRegister(2*norb); c = QuantumCircuit(qr)
c.append(ffsim.qiskit.PrepareHartreeFockJW(norb, nelec), qr)
c.append(ffsim.qiskit.UCJOpSpinBalancedJW(ucj), qr); c.measure_all()
ba = ffsim.qiskit.FfsimSampler(default_shots=200000, norb=norb,
        nelec=nelec, seed=2026).run([c]).result()[0].data.meas

t0 = time.time()
r = diagonalize_fermionic_hamiltonian(
    h, eri, ba, samples_per_batch=spb, norb=norb, nelec=nelec,
    num_batches=1, max_iterations=3,
    sci_solver=partial(solve_sci_batch, spin_sq=0.0, max_cycle=500),
    symmetrize_spin=True, seed=2026)
e = float(r.energy + en); dim = int(np.prod(r.sci_state.amplitudes.shape))
out = {"active_space": "$SPACE", "norb": norb, "nelec": list(nelec),
       "n_qubits": 2*norb, "cas_dim": CAS, "samples_per_batch": spb,
       "energy": e, "subspace_dim": dim, "pct_cas": 100*dim/CAS,
       "wall_s": time.time()-t0}
json.dump(out, open(f"results/sweep/spb{spb:05d}.json", "w"), indent=2)
print(f"spb={spb} E={e:.8f} dim={dim:,} ({100*dim/CAS:.4f}% CAS) {time.time()-t0:.0f}s")
PY

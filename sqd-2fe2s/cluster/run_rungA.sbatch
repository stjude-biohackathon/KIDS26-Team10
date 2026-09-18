#!/bin/bash
#SBATCH --job-name=sqdA-22e16o
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32 --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=sqdA-%j.out --error=sqdA-%j.err
#
# RUNG A  (22e,16o)  32 qubits  19,079,424 determinants
#   Adds the bridging-sulfur 3p superexchange pathway - the thing J depends on.
#   Exact CASCI is still affordable, so SQD gets checked against it at 32 qubits.
#
# Expected: CASCI ladder ~30-40 min, SQD sweeps ~1-2 h.  3 h wall is comfortable.
# This is the rung to run if you only run one.
set -euo pipefail
module load python/3.13 2>/dev/null || true
source "${VENV:-$HOME/venvs/sqd}/bin/activate" 2>/dev/null || true
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}" MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}"
export PYSCF_MAX_MEMORY=$((24 * 1024))
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"; mkdir -p results
T=sto-3g_fe3d+brs3p

echo "=== [A1] exact CASCI ladder (cheapest sector first, checkpointed) ==="
python3 -u src/stage1_classical.py --basis sto-3g --active fe3d+brs3p \
        --ladder-by-sector --tag "$T"

echo "=== [A2] SQD vs exact at 32 qubits ==="
python3 -u src/stage2_sqd.py --tag "$T" --shots 200000 --optimize-ansatz \
        --spb-sweep 200 500 1000 --depol-sweep 0.0 0.01

echo "=== [A3] report ==="
python3 -u src/stage3_report.py --tag "$T" --suffix _opt
echo "RUNG A DONE.  J -> results/stage1_${T}.json .heisenberg.J_cm1"

#!/bin/bash
#SBATCH --job-name=casci-26e18o
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=casci26-%j.out --error=casci26-%j.err
#
# =============================================================================
#  EXACT CASCI at (26e,18o) - the job only the cluster can do
# =============================================================================
#
# 73,410,624 determinants.  The reason this needs a cluster is MEMORY, not cores:
#
#     one CI vector            0.59 GB
#     x max_space=30 Davidson  17.6 GB
#     a 16 GB laptop           cannot hold it
#
# max_space=30 is not optional.  Measured on (10e,10o) against the converged
# energy, the PySCF default max_space=12 lands +0.022 mHa high with one root and
# +0.128 mHa high with three - and the S=0->S=1 gap is only 0.245 mHa.  So a
# smaller Davidson subspace would corrupt J.
#
# What it buys: the exact reference at 36 qubits, which turns the (26e,18o) SQD
# numbers from bare energies into a validated accuracy claim - and gives J at a
# third active space.
#
# Runtime: the (22e,16o) singlet (19.1M dets) took 0.96 h measured.  This is 3.85x
# bigger and the scaling is erratic (measured exponents 1.2-2.5), so honestly:
# somewhere between 5 and 15 hours.  24 h requested.  The ladder is checkpointed
# per sector, cheapest first, so a wall-clock kill still leaves a usable partial
# ladder - three clean spin states are enough to fit J.
# =============================================================================
set -euo pipefail
module load python/3.13 2>/dev/null || true
source "${VENV:-$HOME/venvs/sqd}/bin/activate" 2>/dev/null || true
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}" MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export PYSCF_MAX_MEMORY=$((80 * 1024))
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"; mkdir -p results

python3 -u src/stage1_classical.py --basis sto-3g --active fe3d+brs3p+fe4s \
        --ladder-by-sector --tag sto-3g_fe3d+brs3p+fe4s

echo "DONE. J -> results/stage1_sto-3g_fe3d+brs3p+fe4s.json .heisenberg.J_cm1"
echo "partial ladder (if killed) -> results/stage1_sto-3g_fe3d+brs3p+fe4s_ladder_partial.json"

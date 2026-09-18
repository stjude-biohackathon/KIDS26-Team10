#!/bin/bash
#SBATCH --job-name=sqdC-34e22o
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32 --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=sqdC-%j.out --error=sqdC-%j.err
#
# RUNG C  (34e,22o)  44 qubits  693,479,556 determinants
#   Beyond exact CI: the Davidson would need ~166 GB.  SQD is the only option, which
#   is the regime the method exists for.  No reference of ours, so the energy is
#   reported and compared to published values, not to an exact number we computed.
#
# Expected: stage 1 ~2 min (integrals + CCSD only, CASCI skipped), SQD ~1-1.5 h.
# Holds a 5.5 GB statevector plus ffsim working copies - hence 64 GB.
# Independent of rungs A and B; safe to run in parallel with them.
set -euo pipefail
module load python/3.13 2>/dev/null || true
source "${VENV:-$HOME/venvs/sqd}/bin/activate" 2>/dev/null || true
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}" MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}"
export PYSCF_MAX_MEMORY=$((50 * 1024))
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"; mkdir -p results
T=sto-3g_fe3d+brs3p+2terms3p

echo "=== [C1] integrals + CCSD amplitudes (exact ladder deliberately skipped) ==="
python3 -u src/stage1_classical.py --basis sto-3g --active fe3d+brs3p+2terms3p \
        --skip-casci --tag "$T"

echo "=== [C2] SQD at 44 qubits, no exact reference possible ==="
python3 -u src/stage2_sqd.py --tag "$T" --shots 200000 --optimize-ansatz \
        --spb-sweep 1000 2000 --depol-sweep 0.0 0.01
echo "RUNG C DONE."

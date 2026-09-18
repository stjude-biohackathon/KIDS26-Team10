#!/usr/bin/env bash
# Preflight check — run this on the cluster BEFORE submitting anything.
#
#   bash cluster/preflight.sh
#
# Takes well under a minute and tells you exactly what, if anything, is wrong.
# Everything in this project was developed on macOS / Python 3.13; a cluster is a
# different OS, different BLAS and possibly a different Python, so this script exists
# to surface that difference cheaply instead of after a three-hour queue wait.
#
# Exit code 0 = good to submit.  Non-zero = read the FAIL lines.

set -uo pipefail
cd "$(dirname "$0")/.."

pass=0; warn=0; fail=0
ok   () { echo "  [ OK ]  $*"; pass=$((pass+1)); }
note () { echo "  [WARN]  $*"; warn=$((warn+1)); }
bad  () { echo "  [FAIL]  $*"; fail=$((fail+1)); }

echo "=============================================================="
echo " SQD [2Fe-2S] preflight"
echo " host: $(hostname)   date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "=============================================================="

# ---------------------------------------------------------------- 1. interpreter
echo
echo "1. Python"
PY=${PYTHON:-python3}
if ! command -v "$PY" >/dev/null 2>&1; then
    bad "no '$PY' on PATH. Set PYTHON=/path/to/python3 and re-run."
    echo; echo "SUMMARY: $pass ok, $warn warn, $fail fail"; exit 1
fi
PYV=$("$PY" -c 'import sys;print("%d.%d.%d"%sys.version_info[:3])')
PYMINOR=$("$PY" -c 'import sys;print(sys.version_info[1])')
if [ "$PYMINOR" -ge 11 ]; then ok "Python $PYV at $(command -v "$PY")"
else bad "Python $PYV is too old; need >= 3.11 (developed on 3.13)"; fi

# ---------------------------------------------------------------- 2. packages
echo
echo "2. Packages"
"$PY" - <<'EOF'
import importlib.metadata as md
need = {"pyscf":"2.9","qiskit":"1.3","qiskit-addon-sqd":"0.13","ffsim":"0.0.83",
        "numpy":"1.26","scipy":"1.11","matplotlib":"3.8"}
missing, found = [], {}
for p in need:
    try: found[p] = md.version(p)
    except Exception: missing.append(p)
for p, v in found.items():
    print(f"  [ OK ]  {p} {v}  (need >= {need[p]})")
for p in missing:
    print(f"  [FAIL]  {p} NOT INSTALLED  (need >= {need[p]})")
raise SystemExit(1 if missing else 0)
EOF
if [ $? -eq 0 ]; then pass=$((pass+1)); else bad "install with: $PY -m pip install -r requirements.txt"; fi

# ---------------------------------------------------------------- 3. ffsim patch
echo
echo "3. The ffsim bug our patch works around"
"$PY" - <<'EOF'
import sys, numpy as np
sys.path.insert(0, "src")
try:
    import ffsim.linalg.util as u
except Exception as e:
    print(f"  [FAIL]  cannot import ffsim: {e}"); raise SystemExit(2)
broken = False
try:
    u.unitaries_to_parameters(np.eye(2)[None, :]); print("  [ OK ]  upstream ffsim is already fixed; patch will no-op")
except ValueError:
    broken = True
if broken:
    import ffsim_patch
    if ffsim_patch.PATCHED:
        try:
            u.unitaries_to_parameters(np.eye(2)[None, :])
            print("  [ OK ]  upstream bug present; our patch installed and working")
        except Exception as e:
            print(f"  [FAIL]  patch installed but still broken: {e}"); raise SystemExit(2)
    else:
        print("  [FAIL]  bug present but patch did not install"); raise SystemExit(2)
EOF
if [ $? -eq 0 ]; then pass=$((pass+1)); else bad "ffsim patch problem — optimize=True will fail"; fi

# ---------------------------------------------------------------- 4. resources
echo
echo "4. Resources"
CORES=$("$PY" -c 'import os;print(os.cpu_count())')
ok "cpu_count = $CORES"
MEMKB=$(awk '/MemTotal/{print $2}' /proc/meminfo 2>/dev/null || echo 0)
if [ "$MEMKB" -gt 0 ]; then
    MEMGB=$((MEMKB/1048576))
    echo "  [INFO]  total RAM ${MEMGB} GB"
    [ "$MEMGB" -ge 32 ] && ok "enough RAM for rung A (needs ~8 GB)" \
        || note "rung A wants ~8 GB free; rung B ~24 GB; rung C ~16 GB"
    [ "$MEMGB" -ge 96 ] && ok "enough RAM for rung B (~24 GB Davidson)" \
        || note "rung B may not fit here — it is the one to drop first"
else
    note "cannot read /proc/meminfo (not Linux?); check RAM manually"
fi
for v in OMP_NUM_THREADS MKL_NUM_THREADS; do
    if [ -n "${!v:-}" ]; then ok "$v=${!v}"; else note "$v unset — the sbatch scripts set it"; fi
done

# ---------------------------------------------------------------- 5. smoke test
echo
echo "5. End-to-end smoke test (the real check; ~30-60 s)"
if OUT=$("$PY" - <<'EOF' 2>&1
import sys, time
sys.path.insert(0, "src")
import numpy as np, ffsim_patch                       # noqa: F401
from geometry import build_cluster, to_pyscf_atom
import pyscf.gto, pyscf.scf, pyscf.mcscf, pyscf.ao2mo, pyscf.fci
from pyscf.mcscf import avas
from stage2_sqd import build_circuit, run_sqd, sample_bitstrings

t0 = time.time()
# tiny molecule-free check of the full quantum path: 6 orbitals, 3+3 electrons
rng = np.random.default_rng(0)
norb, nelec = 6, (3, 3)
h = rng.standard_normal((norb, norb)); h = h + h.T
g = rng.standard_normal((norb,) * 4) * 0.1
g = g + g.transpose(1, 0, 2, 3); g = g + g.transpose(0, 1, 3, 2)
g = g + g.transpose(2, 3, 0, 1)
e_exact = float(np.atleast_1d(pyscf.fci.direct_spin1.kernel(h, g, norb, nelec)[0])[0])

t1 = np.zeros((3, 3)); t2 = np.zeros((3, 3, 3, 3)) + 0.05
circ, _ = build_circuit(t1, t2, norb, nelec, 2, True)   # optimize=True -> needs the patch
ba, _ = sample_bitstrings(circ, 2000, norb, nelec, "ffsim", 0.01, 1)
res = run_sqd(h, g, 0.0, ba, norb, nelec, 60, 2, 3, 1, spin_sq=0.0)
print(f"OKSMOKE exact={e_exact:.6f} sqd={res['energy']:.6f} "
      f"dim={res['subspace_dim']} {time.time()-t0:.0f}s")
EOF
); then
    if echo "$OUT" | grep -q OKSMOKE; then
        ok "full path works: $(echo "$OUT" | grep OKSMOKE)"
    else
        bad "smoke test produced no result:"; echo "$OUT" | tail -15
    fi
else
    bad "smoke test crashed:"; echo "$OUT" | tail -20
fi

# ---------------------------------------------------------------- 6. scheduler
echo
echo "6. Scheduler"
if command -v sbatch >/dev/null 2>&1; then
    ok "sbatch found: $(command -v sbatch)"
    command -v sinfo >/dev/null 2>&1 && sinfo -o "  [INFO]  partition %P  nodes %D  cpus %c  mem %m" 2>/dev/null | head -6
else
    note "no sbatch — run the rungs directly with bash instead of submitting"
fi

echo
echo "=============================================================="
echo " SUMMARY: $pass ok, $warn warn, $fail fail"
if [ "$fail" -eq 0 ]; then
    echo " READY. Submit in this order:"
    echo "   sbatch cluster/run_rungA.sbatch     # ~2-2.5 h, do this one first"
    echo "   sbatch cluster/run_rungC.sbatch     # ~1-1.5 h, safe to run in parallel"
    echo "   sbatch cluster/run_rungB.sbatch     # ~3-3.5 h, drop this one if tight"
else
    echo " NOT READY — fix the FAIL lines above. See cluster/README.md for the cut list."
fi
echo "=============================================================="
exit $fail

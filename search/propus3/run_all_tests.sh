#!/bin/sh
# Full validation suite; output of the committed run is validation.log.   sh run_all_tests.sh | tee validation.log
cd "$(dirname "$0")"
PY=${PY:-python3}
T="/usr/bin/time -p"
hdr() { echo; echo "################ $1"; }
hdr "environment"; uname -a; $PY --version; cc --version | head -1; $PY -c "import numpy; print('numpy', numpy.__version__)"; date
hdr "compile p3.c (-Wall -Wextra: no output = no warnings)"
$T cc -O3 -march=native -ffast-math -pthread -Wall -Wextra p3.c -o p3 -lm
hdr "p3struct.selfcheck: class tables and PSD matrix against direct PAF/FFT (v = 239, 133, 219)"
$T $PY -c "
import sys; sys.path.insert(0,'.')
import numpy as np
from p3struct import *
rng = np.random.default_rng(1)
for v, g in ((239,[10]), (133,[26]), (219,[4])):
    st = Structure(v, g); print(v, 'selfcheck', selfcheck(st, rng, 300), '| C', st.C, 'OA', st.A.O, 'OB', st.B.O, 'patterns', [p['rows'] for p in st.feasible_patterns()])
"
hdr "p3test: engine vs independent reference solver on all small structures (v <= 75, pools <= 40000)"
$T $PY p3test.py --vmax 75 --maxpool 40000
hdr "p3validate: repository orders 532 and 876"
$T $PY p3validate.py
hdr "p3cert_test: hit -> matrix -> certificate; repository verify/verify.py as subprocess"
$T $PY p3cert_test.py
hdr "p3plant: planted solution at the full dimensions of v = 239, with negative control"
$T $PY p3plant.py --nb 20000 --nd 300000
hdr "p3scale_test: production configuration (O = 35, C = 17): filters and DFS leaf counts"
$T $PY p3scale_test.py
hdr "p3resume_test: checkpoint / resume"
$T $PY p3resume_test.py
hdr "p3lowe checkc: C energy filter == numpy energy filter"
$T $PY p3lowe.py checkc
hdr "p3lowe study: retention of solution orbits under low-energy cuts"
$T $PY p3lowe.py study
hdr "PSD retention: sampled vs exact"
$T $PY psd_retention_check.py
hdr "Gaussian estimate of the 239 structure (repository estimator, read-only import)"
$T $PY est239.py
hdr "other odd-order multiplier groups of Z_239* (orders 17, 119): exhaustive"
$T $PY p3sweep.py 239 6; $T $PY p3sweep.py 239 2
hdr "v = 67, H of order 3 (Balonin-Djokovic propus structure): exhaustive reduced search"
$T $PY p3sweep.py 67 29
hdr "cost preview (p3plan) for 239"
$T $PY p3plan.py 239 10 --samples 100000
hdr "throughput at the 239 dimensions (single thread, Apple M4): DFS slices (pass 1 only) and the join kernel"
$PY p3problem.py 239 10 --rows 15,1,27 --out "${TMPDIR:-/tmp}/p3_prob239_bench.txt" >/dev/null
for s in 0 21; do $T ./p3 "${TMPDIR:-/tmp}/p3_prob239_bench.txt" -t 8 -countlists -slice $s/64 2>&1 | grep -E "pass 1|real|user"; done
$T ./probebench 40000000
hdr "done"; date

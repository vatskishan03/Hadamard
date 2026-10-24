#!/bin/bash
# quick regression of the whole stack (about 3 minutes, one core); FULL=1 adds the long level-1 test (Z4, r = 3: 7.7e7 instances)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
echo "== t_eq (K-level equations vs G-level PAF)";           python3 t_eq.py
echo "== t_l1 (level 1 vs numpy pair join)";               python3 t_l1.py | grep -E "ALLOK|MISMATCH"
echo "== t_l2 (level 2 vs brute force over all t)";          python3 t_l2.py | grep -v "^[0-9.]*$"
echo "== exhaustive pipeline vs G-level brute force";        for a in "[3] 3 gs" "[2] 3 spence"; do set -- $a; python3 gs_exh.py $1 $2 $3 | tail -1 | cut -c1-200; python3 gbf.py $1 $2 $3 | tail -1; done
echo "== sampler + certificates (gs order 84, spence order 52, bord order 132)"; python3 test_samp.py 2>&1 | grep -E "certificate|result"; python3 test_samp2.py 2>&1 | grep -E "certificate|result"

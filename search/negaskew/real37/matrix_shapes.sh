#!/bin/bash
# the three matrix-H_K shapes of NOTES section 90 (K = Z_p^2, H_K = a Singer subgroup of GL(2, p) of odd order, H = H_K x squares of F_q*): 575 (p = 5, q = 23, order 3), 1183 (p = 13, q = 7, order 21),
# 1083 (p = 19, q = 3, order 45), and the control 275 = Z5^2 x F11 (order 3, the class of skew/1100).  Run in this directory after building qp37r.
# usage: ./matrix_shapes.sh      writes shapes/m*.txt; prints the reduced and the unreduced verdicts
set -e
mkdir -p shapes
python3 qshape37m.py 5 23 '[[2,4],[2,2]]' shapes/m575.txt
python3 qshape37m.py 13 7 '[[10,5],[9,10]]' shapes/m1183.txt
python3 qshape37m.py 19 3 '[[4,6],[3,4]]' shapes/m1083.txt
python3 qshape37m.py 5 11 '[[2,4],[2,2]]' shapes/ctl275.txt
for n in m575 m1183 m1083; do
  echo "== $n reduced";   ${QP:-./qp37r} shapes/$n.txt real=1 l2c=1 first=1 maxprint=1 sym=1 uA=1 tF=1 s3=1 2>&1 | grep -E '^(lists:|STATS|COUNT|SOL)' | cut -c1-200
  echo "== $n unreduced"; ${QP:-./qp37r} shapes/$n.txt real=1 l2c=1 first=1 maxprint=1 2>&1 | grep -E '^(lists:|STATS|COUNT|SOL)' | cut -c1-200
done
echo "== ctl275 (positive control: 9216 normalised solutions unreduced, 24 with every reduction)"
${QP:-./qp37r} shapes/ctl275.txt real=1 l2c=1 maxprint=0 2>/dev/null | grep COUNT
${QP:-./qp37r} shapes/ctl275.txt real=1 l2c=1 maxprint=0 sym=1 uA=1 tF=1 s3=1 2>/dev/null | grep COUNT

#!/bin/bash
# validation of usym / tsym by weighted counts: K = Z9, H_K = <4>, q = 23, 19, 11, with and without the reductions (WCOUNT must be equal).  Run in this directory after building qp37q.
cd "$(dirname "$0")"
mkdir -p shapes
ulimit -v 3500000
export OMP_NUM_THREADS=1
for q in 23 19 11; do
  python3 qshape37q.py '[9]' $q '[4]' shapes/z9_$q.txt >/dev/null
  for opt in "usym=1 tsym=1" "usym=0 tsym=0"; do
    echo "z9 q=$q $opt"
    timeout 300 nice -n 5 ${QP:-./qp37q} shapes/z9_$q.txt l2c=1 maxprint=0 $opt 2>&1 | grep -E "lists built|WCOUNT" | cut -c1-160
  done
done

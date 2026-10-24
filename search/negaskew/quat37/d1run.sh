#!/bin/bash
# defects only on H_K-fixed orbits (d1=1) family; sequential, one core. usage: d1run.sh K q gens tag [extra]      (env TMO: time limit in seconds, default 900; QP: binary, default ./qp37q)
cd "$(dirname "$0")"
K=$1; q=$2; g=$3; tag=$4; shift 4
mkdir -p shapes runs
python3 qshape37q.py "[$K]" $q "[$g]" shapes/d1_$tag.txt > /dev/null
ulimit -v 3500000; export OMP_NUM_THREADS=1
( time timeout ${TMO:-900} nice -n 5 ${QP:-./qp37q} shapes/d1_$tag.txt l2c=1 first=1 maxprint=3 d1=1 usym=1 tsym=1 sym=1 lcap=40000000 "$@" ) > runs/d1_$tag.out 2> runs/d1_$tag.err
echo "EXIT $? $(date -u +%T)" >> runs/d1_$tag.out

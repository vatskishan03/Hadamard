#!/bin/bash
# Sampled runs of qp (thinned lists) on the quaternary plug-in shapes whose lists are too large to enumerate: the UNKNOWN verdicts
# ("sampled without a hit") of NOTES.md section 89.  Every line of a queue is an independent finite job (nice, one core, at most
# 15 or 25 minutes), seeded and so deterministic; a SOL line in an output file would be a hit (certify with plugcert.py).
# usage: ./plug_sampled.sh 1|2|3|299      run in this directory after `gcc -O3 -march=haswell -o qp qp.c -lm`
#   queue 1  plain thinned lists (pL: fraction of the lists of B, C, D; pA = 1: the whole skew list; sym = 1: C <= D)
#   queue 2  block B restricted to the all-E list (bE = 1: no defects; B, C, D are interchangeable and every hit found has such a block)
#   queue 3  bE = 1 with larger samples at the q = 7 shapes (721, 973)
#   299      control: Z_13 x F_23 with <3> (a hit shape), bE = 1, no thinning, C <= D off (sym = 0): 2784 level-1 instances, 2560 solutions, about 1 s
# Logs: runs/sampled<N>.log (one line per job with the STATS line of qp), outputs runs/s<N>_<name>.out / .err.
cd "$(dirname "$0")"
ulimit -v 3500000; export OMP_NUM_THREADS=1
mkdir -p shapes runs
Q=${1:-1}
# shape NAME m q "[H_K generators]"   -> shapes/pNAME.txt
shape() { python3 qshape.py "[$2]" $3 "$4" shapes/p$1.txt > /dev/null; }
shape 737 67 11 "[9]"; shape 869 79 11 "[8]"; shape 943 41 23 "[10]"; shape 779 41 19 "[10]"
shape 759b 33 23 "[4]"; shape 627b 33 19 "[4]"; shape 721 103 7 "[8]"; shape 973 139 7 "[6]"; shape 299 13 23 "[3]"
run() { # run TIMEOUT NAME options...
  tmo=$1; name=$2; shift 2
  echo "$(date -u +%H:%M:%S) start $name $*" >> runs/sampled$Q.log
  timeout $tmo nice -n 5 ./qp shapes/p$name.txt l2c=1 sym=${SYM:-1} maxprint=3 "$@" > runs/s${Q}_$name.out 2> runs/s${Q}_$name.err
  echo "$(date -u +%H:%M:%S) done $name rc=$? $(grep STATS runs/s${Q}_$name.out | cut -c1-260)" >> runs/sampled$Q.log
  if grep -q "^SOL" runs/s${Q}_$name.out; then echo "*** HIT $name" >> runs/sampled$Q.log; fi
}
case $Q in
1) run 900 869 pL=0.05 pA=1 seed=21
   run 900 737 pL=0.06 pA=1 seed=22
   run 900 943 pL=0.10 pA=1 seed=23
   run 900 779 pL=0.04 pA=1 seed=24
   run 900 759b pL=0.01 pA=1 seed=25
   run 900 627b pL=0.01 pA=1 seed=26 ;;
2) run 900 943 bE=1 pL=0.14 pA=1 seed=31
   run 900 779 bE=1 pL=0.026 pA=1 seed=32
   run 900 737 bE=1 pL=0.10 pA=1 seed=33
   run 900 869 bE=1 pL=0.10 pA=1 seed=34
   run 900 759b bE=1 pL=0.007 pA=1 seed=35
   run 900 627b bE=1 pL=0.0065 pA=1 seed=36
   run 900 721 bE=1 pL=0.005 pA=1 seed=37
   run 900 973 bE=1 pL=0.007 pA=1 seed=38 ;;
3) run 1500 721 bE=1 pL=0.03 pA=1 seed=41
   run 1500 973 bE=1 pL=0.03 pA=1 seed=42
   run 1500 721 bE=1 pL=0.06 pA=1 seed=43 ;;
299) SYM=0 run 120 299 bE=1 seed=1 ;;
*) echo "usage: $0 1|2|3|299"; exit 1 ;;
esac

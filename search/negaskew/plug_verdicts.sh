#!/bin/bash
# Exhaustive verdicts of the quaternary plug-in shapes (K = Z_m, q, H_K) of NOTES.md section 89 with the two-level engine qp: the K-level shape, then
# `qp l2c=1` (both levels; level 1 alone gives the number of level-1 instances).  The STATS line has l1_instances (level-1 solutions) and l2_solutions.
# usage: ./plug_verdicts.sh [NAME ...]       run in this directory after `gcc -O3 -march=haswell -o qp qp.c -lm`; no name = every shape except 539
#   NAME is one of 667 851 899 589z31 817a 817b 989 413a 413b 707 539 299 589 (299 and 589 are the hits, 589 on a 1% sample of level 2)
# Expected: l1_instances=0 (EXHAUSTIVELY empty at level 1) for 667 851 899 589z31 817a 817b 989 413a 413b; 707: 290,880 level-1 instances, l2_solutions=0 (about
# 2 minutes); 539: 10,856,160 level-1 instances and l2_solutions=0 (level 1 about 1 minute, level 2 hours: only with the name 539); 299: 36,000 level-1
# instances, 13,056 solutions.  Row-sum verdicts of the same shapes: pal_lists.py (PAL_DRY=1), see README.md.
cd "$(dirname "$0")"
ulimit -v 3500000; export OMP_NUM_THREADS=1
mkdir -p shapes runs
shp() { python3 qshape.py "[$2]" $3 "$4" shapes/v$1.txt > /dev/null; }
verdict() { # NAME m q hK [qp options]
  n=$1; m=$2; q=$3; hk=$4; shift 4
  shp $n $m $q "$hk"
  s=$(date +%s)
  echo "== $n: K = Z_$m, q = $q, H_K = $hk   (qp $*)"
  nice -n 5 ./qp shapes/v$n.txt maxprint=0 "$@" 2>/dev/null | grep -E "^STATS|^COUNT" | cut -c1-250
  echo "   wall $(( $(date +%s)-s )) s"
}
names="$*"; [ -z "$names" ] && names="667 851 899 589z31 817a 817b 989 413a 413b 707 299"
for n in $names; do case $n in
  667) verdict 667 29 23 "[7]" l2c=1 ;;
  851) verdict 851 37 23 "[7]" l2c=1 ;;
  899) verdict 899 29 31 "[7]" l2c=1 ;;
  589z31) verdict 589z31 31 19 "[2]" l2c=1 ;;
  817a) verdict 817a 43 19 "[4]" l2c=1 ;;
  817b) verdict 817b 19 43 "[7]" l2c=1 ;;
  989) verdict 989 43 23 "[4]" l2c=1 ;;
  413a) verdict 413a 7 59 "[]" l2c=1 ;;
  413b) verdict 413b 7 59 "[2]" l2c=1 ;;
  707) verdict 707 101 7 "[5]" l2c=1 ;;
  539) verdict 539 49 11 "[2]" l2c=1 ;;
  299) verdict 299 13 23 "[3]" l2c=1 ;;
  589) verdict 589 19 31 "[7]" l2c=1 l2p=0.01 first=1 seed=3 ;;
  *) echo "unknown name $n" ;;
esac; done

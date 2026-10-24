#!/bin/bash
# restartable production worker for the real GS-type plug-in sampler.  One chunk = one level-1 + level-2 round = one python process (a restart loses at most the round in progress).
# usage: gs_worker.sh NAME MODE KMODS R HG T ALLOCS LA LC SEEDBASE
#   e.g. gs_worker.sh w1 gs '[15]' 7 '[]' null "2:2:2:2@300,3e5,..."  (ALLOCS: comma list of  n0[.m0]:n1[.m1]:n2[.m2]:n3[.m3][/R0:R1:R2:R3][@LA,LC]  -- see gs_samp.py; '@' list sizes optional)
# GS_SECS (default 900) = seconds per python process; a restart loses at most that much.
# stops when runs/NAME.STOP exists or a HIT file runs/NAME*.HIT_* exists.  Uses ONE core: run several workers with different NAMEs / SEEDBASE only if cores are free.
cd "$(dirname "$0")"
mkdir -p runs
export OMP_NUM_THREADS=1
ulimit -v 3500000
NAME=$1; MODE=$2; KM=$3; R=$4; HG=$5; T=$6; ALLOCS=$7; LA=$8; LC=$9; SB=${10}
SECS=${GS_SECS:-900}      # seconds per python process (several rounds each; patterns are computed once per process)
IFS=';' read -ra AL <<< "$ALLOCS"      # allocations separated by ';' (the '@LA,LC' list sizes contain commas)
i=0
while [ ! -e runs/$NAME.STOP ] && ! ls runs/${NAME}.HIT_* >/dev/null 2>&1; do
  a=${AL[$((i % ${#AL[@]}))]}
  nice -n 5 python3 gs_samp.py "$MODE" "$KM" $R "$HG" "$T" "$a" $LA $LC $SECS $((SB + i)) runs/$NAME >> runs/$NAME.log 2>&1
  i=$((i+1))
done
echo "worker $NAME finished after $i chunks at $(date -u)" >> runs/$NAME.log

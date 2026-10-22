#!/bin/bash
# restartable production worker for a plain plug-in shape (cyclic K).  One chunk = one level-1/level-2 round = one python process.
# usage: prod_worker.sh NAME KMODS Q HG ALLOCS LC SEEDBASE      (e.g. prod_worker.sh p418w1 '[19]' 11 '[]' "1:1:2:3@1.2e7,1:1:2:3@1.2e7,1:1:1:4@1e7" 12000000 $(date +%s))
# stops when runs/NAME.STOP exists or any runs/NAME*.FOUND / runs/<shape>.FOUND exists
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
ulimit -v 3500000
NAME=$1; KM=$2; Q=$3; HG=$4; ALLOCS=$5; LC=$6; SB=$7; SHAPE=${8:-$1}
IFS=',' read -ra AL <<< "$ALLOCS"
i=0
while [ ! -e runs/$NAME.STOP ] && ! ls runs/${SHAPE}*.FOUND >/dev/null 2>&1; do
  a=${AL[$((i % ${#AL[@]}))]}
  nice -n 5 python3 samp_l1.py "$KM" $Q "$HG" "$a" $LC 1 $((SB + i)) runs/$NAME >> runs/$NAME.log 2>&1
  i=$((i+1))
done
echo "worker $NAME finished after $i chunks at $(date -u)" >> runs/$NAME.log

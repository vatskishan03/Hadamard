#!/bin/bash
# runjoin.sh INPUT NS FIRST LAST STEP OUTDIR [MAXPAIRS=6e7] [PAR=0] [NPAR=1]
# Restart-safe driver of rj: bucket ranges [i, i+STEP) are run one after another, each leaving OUTDIR/b<i>_<j>.{out,err,done};
# only ranges with (range index mod NPAR) == PAR are run (so NPAR=2 instances with PAR=0,1 share the work on 2 cores).
# Solutions are appended to OUTDIR/sol_all.txt; the driver keeps going unless OUTDIR/STOP exists.
set -u
IN=$1; NS=$2; FIRST=$3; LAST=$4; STEP=$5; OUT=$6; MP=${7:-6e7}; PAR=${8:-0}; NPAR=${9:-1}
mkdir -p "$OUT"
ulimit -v 3500000
export OMP_NUM_THREADS=1
HERE="$(cd "$(dirname "$0")" && pwd)"
i=$FIRST; c=0
while [ "$i" -lt "$LAST" ]; do
  j=$((i+STEP)); [ "$j" -gt "$LAST" ] && j=$LAST
  tag="b${i}_${j}"
  if [ $((c % NPAR)) -eq "$PAR" ] && [ ! -f "$OUT/$tag.done" ] && [ ! -f "$OUT/STOP" ]; then
    nice "$HERE/rj" "$IN" "$NS" "$LAST" "$MP" "$i" "$j" > "$OUT/$tag.out" 2> "$OUT/$tag.err" || { echo "rj failed at $tag"; exit 1; }
    touch "$OUT/$tag.done"
    if grep -q '^SOL' "$OUT/$tag.out"; then grep '^SOL' "$OUT/$tag.out" >> "$OUT/sol_all.txt"; echo "SOLUTION in $tag"; fi
  fi
  i=$j; c=$((c+1))
done

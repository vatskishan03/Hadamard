#!/bin/bash
# resumable driver for tri2 (search/NOTES.md section 92; README.md of this directory): usage  ./run_tri2.sh TAG v k g mode red M r0 step cnt maxsol seed
# appends to runs/TAG.out (SOL and STAT lines) and runs/TAG.err; resumes at the number of STAT lines already in runs/TAG.out (partition index idx0 of tri2).
# One core (the caller chooses it, e.g. with taskset), nice, 3.5 GB; detach long runs with  setsid nohup ./run_tri2.sh ... &  and record the PID (echo $! > runs/TAG.pid).
cd "$(dirname "$0")" || exit 1
mkdir -p runs
TAG=$1; shift
OUT=runs/$TAG.out
DONE=0; [ -f $OUT ] && DONE=$(grep -c '^STAT' $OUT)
ulimit -v 3500000
echo "$(date -u +%FT%TZ) start/resume $TAG at partition $DONE" >> runs/$TAG.log
exec nice -n 5 env OMP_NUM_THREADS=1 ./tri2 "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "$DONE" "${10}" "${11}" >> $OUT 2>> runs/$TAG.err

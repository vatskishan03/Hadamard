#!/bin/bash
# usage: run_exh.sh NAME SECS 'K' q 'hunits' SLOTFILE TOPN ORDER ENGINE   (env EX_* passed through)   -- resumable (runs/NAME.slots.jsonl, runs/NAME.seen); pid in runs/NAME.pid
cd "$(dirname "$0")/.."
mkdir -p runs
ulimit -v 3500000
export OMP_NUM_THREADS=1
echo $$ > runs/$1.pid
exec nice -n 5 python3 -u exh38.py "$@"

#!/bin/bash
# usage: run_job.sh NAME SEED SECS 'K' q 'hunits' TOPN LB LC LD [ENGINE]      (the slot file est/NAME.jsonl must exist)
# runs prod_dedup.py in search/skewplug (resumable through runs/NAME.jsonl and runs/NAME.seen)
cd "$(dirname "$0")/.."
mkdir -p runs
ulimit -v 3500000
export OMP_NUM_THREADS=1 PROD_EREL=5 PROD_MAXOUT=${PROD_MAXOUT:-20000} PROD_L2CAP=${PROD_L2CAP:-3000}
exec nice -n 5 python3 -u prod_dedup.py "$1" "$2" "$3" "$4" "$5" "$6" EST "$7" "$8" "$9" "${10}" "${11:-skf}"

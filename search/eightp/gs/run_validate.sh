#!/bin/bash
# sequential validation battery (one core): gs, spence, bord;  resume-safe
cd "$(dirname "$0")"; mkdir -p runs
export OMP_NUM_THREADS=1
ulimit -v 3500000
for m in gs spence bord; do
  nice -n 5 python3 -u validate.py $m runs/val_$m.jsonl 60 >> runs/val_$m.log 2>&1
done
echo done >> runs/val_done.txt

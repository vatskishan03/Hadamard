#!/bin/bash
# second validation battery with the final binaries: 22 shapes per mode, >= 8 nonzero.  Resume-safe.
cd "$(dirname "$0")"; mkdir -p runs
export OMP_NUM_THREADS=1
ulimit -v 3500000
for m in gs spence bord; do
  nice -n 8 python3 -u validate.py $m runs/val2_$m.jsonl 40 22 8 >> runs/val2_$m.log 2>&1
done
echo done >> runs/val2_done.txt

#!/bin/bash
# exhaustive d = 3 spin structures (A, B, muB, mu^2B) over Z_v, v = 2p-1, all pairs H < K of index 3, O_H <= OMAX; one core; resumable per v (.done)
cd "$(dirname "$0")"
mkdir -p runs
export OMP_NUM_THREADS=1
ulimit -v 3000000
OMAX=${OMAX:-26}
for v in 1197 1477 1677 2061; do
  if [ ! -e runs/v$v.O$OMAX.done ]; then
    nice -n 5 python3 spin_drv.py $v $OMAX 2 > runs/v$v.O$OMAX.log 2>&1 && touch runs/v$v.O$OMAX.done
  fi
done
echo finished > runs/run_spin.finished

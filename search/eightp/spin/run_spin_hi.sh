#!/bin/bash
# d = 3 spin structures (A, B, mu B, mu^2 B) for OMIN <= O_H <= OMAX (continuation of run_spin.sh), one core, resumable per v
cd "$(dirname "$0")"
mkdir -p runs
export OMP_NUM_THREADS=1
ulimit -v 3000000
OMAX=${OMAX:-30}
OMIN=${OMIN:-27}
for v in 1197 1477 2061 1677; do
  if [ ! -e runs/v$v.O$OMIN-$OMAX.done ]; then
    nice -n 5 python3 spin_drv.py $v $OMAX $OMIN > runs/v$v.O$OMIN-$OMAX.log 2>&1 && touch runs/v$v.O$OMIN-$OMAX.done
  fi
done
echo finished > runs/run_spin_hi.$OMIN-$OMAX.finished

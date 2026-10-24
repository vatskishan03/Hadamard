#!/bin/bash
# exhaustive d = 2 and d = 4 block-symmetric Spence structures over Z_v (spinD_drv.py), all (H, K) with O_H <= OMAX; one core; resumable per (v, d) (.done)
cd "$(dirname "$0")"
mkdir -p runs
export OMP_NUM_THREADS=1
ulimit -v 3000000
OMAX=${OMAX:-26}
OMIN=${OMIN:-2}
for d in 4 2; do
for v in 1197 1437 1477 1677 2061 2173; do
  if [ ! -e runs/D$d.v$v.O$OMIN-$OMAX.done ]; then
    nice -n 5 python3 spinD_drv.py $v $d $OMAX $OMIN > runs/D$d.v$v.O$OMIN-$OMAX.log 2>&1 && touch runs/D$d.v$v.O$OMIN-$OMAX.done
  fi
done
done
echo finished > runs/run_spinD.$OMIN-$OMAX.finished

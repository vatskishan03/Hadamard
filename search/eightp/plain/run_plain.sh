#!/bin/bash
# exhaustive plain cyclic Spence shapes (odd v = 2p-1), all subgroups H <= Z_v^* with 2..OMAX orbits (solver plain_spence3.py); resumable per v (a finished v writes runs/vV.O<OMAX>.done)
cd "$(dirname "$0")"
mkdir -p runs
export OMP_NUM_THREADS=1
ulimit -v 3500000
OMAX=${OMAX:-22}
for v in 1437 1477 2173 1677 2061 1197; do
  if [ ! -e runs/v$v.O$OMAX.done ]; then
    nice -n 5 python3 plain_spence3.py $v $OMAX 2 > runs/v$v.O$OMAX.log 2>&1 && touch runs/v$v.O$OMAX.done
  fi
done
echo finished > runs/run_plain.finished

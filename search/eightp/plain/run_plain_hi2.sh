#!/bin/bash
# one-core queue: v = 2061, 2173 at O = 23..24 with the final plain_spence3.py (exact additive residues mod nparts; chunked pairs); an earlier solver ran out of memory (MemoryError) at these two v.  WAITPID = pid to wait for.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
ulimit -v 3500000
mkdir -p runs
while kill -0 "${WAITPID:-0}" 2>/dev/null; do sleep 20; done
for v in 2061 2173; do
  [ -e runs/r4_${v}_23_24.done ] || (nice -n 5 python3 plain_spence3.py $v 24 23 1e5 > runs/r4_${v}_23_24.log 2>&1 && touch runs/r4_${v}_23_24.done)
done
echo finished > runs/run_plain_hi2.finished

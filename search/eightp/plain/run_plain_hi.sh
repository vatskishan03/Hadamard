#!/bin/bash
# one-core queue: O = 23..24 for v = 1197, 1677, 2061, 2173 and the re-verification of v = 1477 (O <= 22) and of the non-cyclic Z_3^2 x Z_133 sample with the corrected solver plain_spence3.py (search/NOTES.md section 92).
# The first solver split large pair buckets by (hash mod nparts) with nparts not a power of two; the int64 hash sums wrap around, so the residues were inconsistent => a false EMPTY was possible for buckets > 4e7 pairs.
# Resumable via .done files; logs runs/r3_*.log (solver times of the recorded logs: v = 1477 50 s, v = 1197 218 s, v = 1677 412 s)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
ulimit -v 3500000
mkdir -p runs
job() { # tag args...
  tag=$1; shift
  [ -e runs/r3_$tag.done ] || (nice -n 5 python3 "$@" > runs/r3_$tag.log 2>&1 && touch runs/r3_$tag.done)
}
job 1477_2_22  plain_spence3.py 1477 22 2 1e5
job 1197_23_24 plain_spence3.py 1197 24 23 1e5
job 1677_23_24 plain_spence3.py 1677 24 23 1e5
job 2061_23_24 plain_spence3.py 2061 24 23 1e5
job 2173_23_24 plain_spence3.py 2173 24 23 1e5
job nc_z3z3_133 plain_nc3.py z3z3 133 40000 24 1e5
echo finished > runs/run_plain_hi.finished

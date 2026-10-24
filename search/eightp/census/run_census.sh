#!/bin/bash
# the census of the real plug-in shapes for n = 2p, p = 599, 719, 739, 839, 1031, 1087 (search/NOTES.md section 92): shapes G = K x F_r with 2..40 orbits and the exact identity / row-sum filters of shapes4.dp_filter, for
#   cens_<n>   mode b (Spence, |G| = n - 1), r >= 5, cyclic multiplier groups     censP_<n>  the same with product groups
#   cens3_<n>  r >= 3, cyclic                                                      cens3P_<n> r >= 3, product groups
#   censC_<n>  mode c (bordered, Q = Z_2, 4, 5, 6, 7 ...)
# and a comparison of the outputs (runs/) with the recorded ones (../data/census).  About 100 s on one core.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
mkdir -p runs
for n in 1198 1438 1478 1678 2062 2174; do
  python3 census_a.py $n runs/cens_$n.jsonl 2 40
  python3 census_a.py $n runs/censP_$n.jsonl 2 40 product
  python3 census_a.py $n runs/cens3_$n.jsonl 2 40 cyclic z3,z2z2 b 3
  python3 census_a.py $n runs/cens3P_$n.jsonl 2 40 product z3,z2z2 b 3
  python3 census_a.py $n runs/censC_$n.jsonl 2 40 cyclic z2,z4,z5,z6,z7 c
done
same=0; diff=0
for f in runs/cens*_*.jsonl; do b=$(basename $f); if cmp -s $f ../data/census/$b; then same=$((same+1)); else diff=$((diff+1)); echo "DIFFERENT $b"; fi; done
echo "census: $same files identical to ../data/census, $diff different"

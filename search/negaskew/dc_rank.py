"""Row-sum test (dcnp.plain_dc) of every ranked shape with n = 15..21 orbits and estimate >= -9 bits.
usage: python3 enum_shapes.py 21 2499 && python3 rank_shapes.py && python3 dc_rank.py     (reads runs/rank_n15_21.json)
Prints (estimate, v, Table 5 entry m, n, |H|, generators, 'DC' True/False, orbit sizes).  Only 101 (n = 21), 299 (n = 15),
589 <568> and <586> (n = 17), 821 (n = 21), 473 (n = 21) and 817 (n = 17) have a row-sum pattern; 337, 461, 401, 487,
523, 491, 661, 659 and 739 (all Table 5 entries m >= 4) die by it: an orbit size s >= 21 gives R = u0 + s w and no
combination has sum |R|^2 = 4v (NOTES.md section 88)."""
import json
from qdrv import shape
from dcnp import plain_dc
r=json.load(open('runs/rank_n15_21.json'))
for s,v,m,n,o,g,dc in r:
    if s< -9: continue
    sh=shape(v,g)
    d=plain_dc(v,sh['sz'],sh['neg'],n)
    print(s,v,m,n,o,g,'DC',d,'sizes',sorted(set(sh['sz'])),flush=True)

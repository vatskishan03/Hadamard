"""Ranks the multiplier shapes with 15 to 21 orbits (beyond the four-list join of qskew) by the Gaussian estimate of est.py,
and tests the row-sum condition (dcnp.plain_dc) of those with estimate >= 5 bits and v <= 1500.
usage: python3 enum_shapes.py 21 2499 && python3 rank_shapes.py      writes runs/rank_n15_21.json, prints the 70 best
Each row is (estimate in bits, v, Table 5 entry m, orbits n, |H|, generators, row-sum pattern exists).  The shapes with
n <= 13 are the sweep of sweep.py and are left out.  Of the shapes found with an estimate >= 5 bits only v = 101 (n = 21,
<95> of order 5) has a row-sum pattern; 337, 461, 401 and 487 die by the row sums (NOTES.md section 88)."""
import json,os
from est import est
from qdrv import shape
from dcnp import plain_dc
d=json.load(open('enum_n21_v2499.json'))
out=[]
for v,(m,sh) in sorted(d.items(),key=lambda x:int(x[0])):
    v=int(v)
    for n,o,g in sh:
        if n<=13: continue   # n <= 13 is the sweep of sweep.py
        b,l,s=est(v,g)
        dc=None
        if s>=5 and v<=1500:
            sh_=shape(v,g); dc=plain_dc(v,sh_['sz'],sh_['neg'],n)
        out.append((round(s,1),v,m,n,o,g,dc))
out.sort(key=lambda r:-r[0])
os.makedirs('runs',exist_ok=True)
json.dump(out,open('runs/rank_n15_21.json','w'))
print(len(out))
for r in out[:70]: print(r)

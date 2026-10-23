"""Number of B that pass the row-sum test (passDC) for each DC-feasible spin shape with n = 21 orbits, from the exact
distribution of the row sums (costest.dpB / dpA); qspin_fast spends about 0.2 us per such B plus a setup of 45 to 80 s.
usage: python3 spin_est.py && python3 spin_filter.py && python3 spin_cost.py      (reads runs/spin_rank_dc.json, n = 21 only)
Example: v = 589, H = <574>, mu = 5: passDC 6.15e7 (the measured value is 61,471,872)."""
import json, sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..'))
from qdrv import shape
from costest import dpB, dpA, norm_hist
rows=[r for r in json.load(open('runs/spin_rank_dc.json')) if r[1]==21]
for sp,n,v,m,o,g,mu,K in rows:
    d=shape(v,g); cB=dpB(v,d['sz'],d['n']); cA=dpA(v,d['sz'],d['neg'],d['n'])
    hA=norm_hist(cA,v); tot=4*v
    S=2*v+1; xs,ys=np.nonzero(cB); 
    pass_dc=0.0
    for x,y in zip(xs,ys):
        nm=3*((x-v)**2+(y-v)**2)
        if nm<=tot and (tot-nm) in hA: pass_dc+=cB[x,y]
    print(v,o,g,'mu',mu,'est',sp,'passDC %.3g'%pass_dc,flush=True)

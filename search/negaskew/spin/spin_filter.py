"""Row-sum filter of the ranked spin shapes: keeps those with a pattern (|R_A|^2 + 3 |R_B|^2 = 4v, dcnp.spin_dc).
usage: python3 spin_est.py && python3 spin_filter.py     reads runs/spin_rank.json, writes runs/spin_rank_dc.json
At n = 17 only 2 of the 38 shapes pass (v = 817, two groups); at n = 19 none of 13 (NOTES.md section 88)."""
import json,time,os,sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..'))
from qdrv import shape
from dcnp import spin_dc
rk=json.load(open('runs/spin_rank.json'))
keep=[]
for sp,n,v,m,o,g,mu,K in rk:
    t0=time.time()
    d=shape(v,g)
    ok=spin_dc(v,d['sz'],d['neg'],d['n'])
    print(v,o,g,mu,'est',sp,'DC',ok,round(time.time()-t0,1),flush=True)
    if ok: keep.append((sp,n,v,m,o,g,mu,K))
json.dump(keep,open('runs/spin_rank_dc.json','w'))
print(len(keep),'of',len(rk))

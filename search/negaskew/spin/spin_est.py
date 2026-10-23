"""Ranks the spin shapes with n >= 14 orbits (runs/spin_shapes_n<NMAX>.json of enum_spin.py, NMAX = 17 or 21) by the spin estimate
sp = (3 (n-1) + 1) + l/3 bits, where l is the Gaussian equation probability of est.py; the shapes already in runs/spin_sweep.jsonl
are left out.  usage: python3 spin_est.py [NMAX]   -> runs/spin_rank.json (read by spin_filter.py), prints the ranking (NMAX default 17).
The estimate over-predicts by 2 to 6 bits when calibrated on v = 7 and 13 (12.5 / 22.3 against 10.8 / 16.8 bits exact),
and by at least 9 bits at the five exhaustive n = 21 shapes (NOTES.md section 89)."""
import json,os,sys
nmax=int(sys.argv[1]) if len(sys.argv)>1 else 17
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..'))
from est import est
sh=json.load(open(f'runs/spin_shapes_n{nmax}.json'))
done=set()
if os.path.exists('runs/spin_sweep.jsonl'):
    for l in open('runs/spin_sweep.jsonl'):
        r=json.loads(l); done.add((r['v'],tuple(r['gens']),r['mu']))
out=[]
for n,v,m,o,g,mu,K in sh:
    if n<14 or (v,tuple(g),mu) in done: continue
    b,l,s=est(v,g)
    sp=(3*(n-1)+1)+l/3.0
    out.append((round(sp,1),n,v,m,o,g,mu,K))
out.sort(reverse=True)
json.dump(out,open('runs/spin_rank.json','w'))
for r in out: print(r)

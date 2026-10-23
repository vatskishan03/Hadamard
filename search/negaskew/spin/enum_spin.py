"""Lists the spin shapes (H, mu): for each census-undecided v (the list of enum_shapes.py) and each multiplier group H of enum_n21_v2499.json
with at most NMAX orbits, the classes of mu in Z_v* with mu not in H and mu^3 in H (C = B o mu, D = B o mu^2).
usage: (in ..)  python3 enum_shapes.py 21 2499        then here:  python3 enum_spin.py 17      -> runs/spin_shapes_n17.json
Each row is (n orbits, v, Table 5 entry m, |H|, generators of H, mu, order of <H, mu>)."""
import json,sys,os
from math import gcd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..'))
from enum_shapes import gen_sub, norb
nmax=int(sys.argv[1])
d=json.load(open(os.path.join(HERE,'..','enum_n21_v2499.json')))
out=[]
for v,(m,sh) in d.items():
    v=int(v)
    units=[h for h in range(2,v) if gcd(h,v)==1]
    for n,o,g in sh:
        if n>nmax: continue
        H=gen_sub(v,g)
        seen=set()
        for mu in units:
            if mu in H: continue
            if pow(mu,3,v) not in H: continue
            K=gen_sub(v,list(g)+[mu])
            if K in seen: continue
            seen.add(K)
            out.append((n,v,m,o,g,mu,len(K)))
out.sort()
os.makedirs('runs',exist_ok=True)
json.dump(out,open(f'runs/spin_shapes_n{nmax}.json','w'))
from collections import Counter
print(len(out),Counter(x[0] for x in out))

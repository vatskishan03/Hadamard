"""Paley shapes among the multiplier shapes with at most 21 orbits: v = K q with q = 3 mod 4 prime >= 7, gcd(K, q) = 1, and H
containing {1} x QR(q); kappa = number of K-orbits of the image H_K.
usage: python3 enum_shapes.py 21 2499 && python3 paley_shapes.py      writes runs/paley_shapes.json (read by idfeas.py)
Rows: (v, Table 5 entry m, q, K, n orbits, |H|, generators, kappa), sorted by kappa."""
import json, sys, os
from math import gcd
from sympy import factorint, isprime
from enum_shapes import gen_sub, norb
d=json.load(open('enum_n21_v2499.json'))
t5=json.load(open('t5.json'))
out=[]
for v,(m,sh) in d.items():
    v=int(v)
    f=factorint(v)
    for p,e in f.items():
        if e!=1 or p%4!=3 or p<7: continue
        K=v//p
        if gcd(K,p)!=1: continue
        # QR subgroup of Z_v^*: elements = 1 mod K, square mod p
        # generator: crt
        for n,o,g in sh:
            H=gen_sub(v,g)
            # H contains {1}x QR(p)?  elements h with h%K==1 and h%p square
            qr=[h for h in range(1,v) if h%K==1 and pow(h%p,(p-1)//2,p)==1 and gcd(h,v)==1]
            if all(h in H for h in qr):
                # H_K = image of H mod K ; kappa = number of K-orbits under H mod K
                HK=set(h%K for h in H)
                kappa=norb(K,HK) if K>1 else 1
                out.append((v,m,p,K,n,o,g,kappa))
out.sort(key=lambda r:(r[7],r[0]))
for r in out: print(r)
os.makedirs('runs',exist_ok=True)
json.dump(out,open('runs/paley_shapes.json','w'))
print(len(out))

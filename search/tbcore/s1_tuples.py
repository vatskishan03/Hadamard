"""Step 1 of the n = 955 search: the exact join at the level of Z29 (the five Z29-level class values of each block).

G = Z2^3 x Z29, H = <(tau, 16)>, Q = Z2^3, T = the seven nontrivial characters (pattern 2: row sums 21, 21, 6, 6; compressions
yA = yB = (0, 3), yC = yD = (-1, 1) on Q). Each block is summed over the Z2^3 coordinates (z-compression, 5 H-orbit values on
Z29) and over Z29 (Q-compression); the A- and C-type joint classes are listed with their multiplicities, the PAF of the Z29
compression is reduced to 3 numbers and the four lists are joined exactly. The tuples are written to OUT/tuples.json,
sorted by decreasing weight (the index of a tuple in this list is the index used by s1_prod.py).

usage: python3 s1_tuples.py [--out DIR]
"""
import argparse, itertools, math, time
from math import comb
from collections import defaultdict
import numpy as np
import tbio

ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("--out", default=None, help="output directory (default: current directory)")
args = ap.parse_args()
out = tbio.outdir(args.out)

p=29; u=16; g=2
# orbits of <u> on Z29*: O_t = {g^(t+4i)}
ind = {pow(g,k,p):k for k in range(28)}
O = [[pow(g,t+4*i,p) for i in range(7)] for t in range(4)]
assert all(pow(u,7,p)==1 for _ in [0])
# y-vector on Z29 from Z29-level class vector (y0, y_t)
def yvec(y):
    v=np.zeros(p,dtype=np.int64); v[0]=y[0]
    for t in range(4):
        for z in O[t]: v[z]=y[t+1]
    return v
def paf(v):
    return np.array([np.dot(v,np.roll(v,-w)) for w in range(p)])
Svals=[-7,-5,-3,-1,1,3,5,7]; mult={S:comb(7,(7+S)//2) for S in Svals}
# C-type
Ccls=defaultdict(int)   # (Zvec, Qvec) -> count
for c0 in (1,-1):
  for cE in (1,-1):
    for ct in itertools.product((1,-1),repeat=4):
      for S in itertools.product(Svals,repeat=4):
        z=(c0+7*cE,)+tuple(ct[t]+S[t] for t in range(4))
        q=(c0+7*sum(ct), cE+sum(S))
        m=1
        for t in range(4): m*=mult[S[t]]
        Ccls[(z,q)]+=m
Acls=defaultdict(int)
for a in (1,-1):
  for al in itertools.product((1,-1),repeat=2):
    for S0 in Svals:
      for S1 in Svals:
        z=(7*a, al[0]+S0, al[1]+S1, al[0]+S0, al[1]+S1)
        q=(14*(al[0]+al[1]), a+2*S0+2*S1)
        Acls[(z,q)]+=mult[S0]*mult[S1]
print("C joint classes",len(Ccls),"A joint classes",len(Acls))
# restrict to Q tuple of pattern 2
QA=(0,3); QC=(-1,1)
A=[(z,c) for (z,q),c in Acls.items() if q==QA]
C=[(z,c) for (z,q),c in Ccls.items() if q==QC]
print("A classes with Q=(0,3):",len(A),"total",sum(c for z,c in A)," C classes:",len(C),"total",sum(c for z,c in C))
# join at Z29 level: sum paf = 0 for w != 0 ; total sum squares 954
pa={z:paf(yvec(z)) for z,c in A}
pc={z:paf(yvec(z)) for z,c in C}
# reduce PAF vector to the 2 classes (w in O0 U O2 , w in O1 U O3) plus w=0
def red(v): return (int(v[0]), int(v[O[0][0]]), int(v[O[1][0]]))
for d in (pa,pc):
    for z in d: d[z]=red(d[z])
# A+B lists
AB=defaultdict(list)
for za,ca in A:
    for zb,cb in A:
        key=tuple(x+y for x,y in zip(pa[za],pa[zb]))
        AB[key].append((za,zb,ca*cb))
tot=0; tuples=[]
CD=defaultdict(list)
for zc,cc in C:
    for zd,cd in C:
        key=tuple(x+y for x,y in zip(pc[zc],pc[zd]))
        CD[key].append((zc,zd,cc*cd))
tgt=(954,0,0)
for k1,l1 in AB.items():
    k2=tuple(t-x for t,x in zip(tgt,k1))
    if k2 in CD:
        for (za,zb,w1) in l1:
            for (zc,zd,w2) in CD[k2]:
                tuples.append((za,zb,zc,zd,w1*w2))
print("joint feasible tuples:",len(tuples),"mass 10^%.2f"%math.log10(sum(t[4] for t in tuples)))
tuples.sort(key=lambda t:-t[4])
tbio.save_tuples(out, tuples)
print("wrote", tbio.tuples_path(out), " sha256 of the canonical list", tbio.tuples_digest(tuples))
for t in tuples[:8]: print(t[:4], "%.3g"%t[4])

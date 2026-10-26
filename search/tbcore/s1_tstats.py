"""Step 2: for every tuple of OUT/tuples.json the sizes of its four lists.

For tuple (zA, zB, zC, zD, w): |A| = |B| is the number of A-type sequences (exhaustive enumeration of the 2^17) with Q-compression
(0, 3) and Z29-class zA (resp. zB); the 'fibre' of C (resp. D) is the number of C-type sequences with Q-compression (-1, 1) and
Z29-class zC (zD), counted exactly. Writes OUT/trows.npy, rows (index, |A|, |B|, fibre C, fibre D, weight).

usage: python3 s1_tstats.py [--out DIR]
"""
import argparse
import numpy as np
from collections import Counter
import tbio
from s1_lists import *
from s1_prepare import C_combos

ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("--out", default=None, help="output directory (default: current directory)")
args = ap.parse_args()
out = tbio.outdir(args.out)

tuples=tbio.load_tuples(out)
f0,f1=enum_A(); XA=build(f0,f1); z5,q0,q1=zq_compress(XA)
m=(q0==0)&(q1==3)
cntA=Counter(map(tuple,z5[m].tolist()))
cC={}
def fibreC(z):
    if z not in cC: cC[z]=sum(c[4] for c in C_combos(z,(-1,1)))
    return cC[z]
rows=[]
for ti,(za,zb,zc,zd,w) in enumerate(tuples):
    rows.append((ti,cntA[za],cntA[zb],fibreC(zc),fibreC(zd),w))
rows=np.array(rows,dtype=float)
print("distinct |A| values:",sorted(set(int(x) for x in rows[:,1])))
print("fibre C range: %.3g .. %.3g"%(rows[:,3].min(),rows[:,3].max()))
AB=rows[:,1]*rows[:,2]
order=np.argsort(-AB)
print("top by |A||B|:")
for i in order[:6]: print(rows[i][:6])
print("median |C|: %.3g"%np.median(rows[:,3]))
# number of tuples by |A||B|
print(Counter((int(a),int(b)) for a,b in rows[:,1:3]).most_common(8))
np.save(tbio.trows_path(out),rows)
print("wrote", tbio.trows_path(out))

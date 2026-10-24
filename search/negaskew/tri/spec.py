"""coset-DFT supports and magnitudes of the blocks of a solution file (blocks are functions on the e = (v-1)/k cosets): python3 spec.py v "[H generators]" FILE [MAXSHOW]
(v = 1 + k(k+1)/2 prime, H the subgroup of order k of Z_v*; FILE: SOL lines of ../qskew or tri_join; the blocks are quaternary exponent strings, one digit per H-orbit)"""
import collections, numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))      # ../qdrv.py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qdrv import shape
from sympy import primitive_root
v=int(sys.argv[1]); gens=eval(sys.argv[2]); fn=sys.argv[3]; maxshow=int(sys.argv[4]) if len(sys.argv)>4 else 5
d=shape(v,gens); n=d['n']; H=set(d['H']); k=len(H); e=(v-1)//k
g=primitive_root(v)
# coset index of each orbit (orbits 1..n-1)
cos=[None]*n
for o in range(1,n):
    r=d['reps'][o]
    for j in range(e):
        if (pow(g,j,v)*1) % v in [(r*pow(h,-1,v))%v for h in H]: cos[o]=j; break
I=np.array([1,1j,-1,-1j])
om=np.exp(2j*np.pi/e)
supp=collections.Counter(); shown=0; mags=collections.Counter()
for l in open(fn):
    if not l.startswith('SOL'): continue
    E=[[int(c) for c in p] for p in l.split()[1:]]
    row=[]
    for X,ex in enumerate(E):
        c=np.zeros(e,complex)
        for o in range(1,n): c[cos[o]]=I[ex[o]]
        ch=np.array([c @ om**(-m*np.arange(e)) for m in range(e)])/e   # chat(m)
        S=c.sum()
        nz=tuple(m for m in range(e) if abs(ch[m])>1e-9)
        row.append((X,ex[0],np.round(S,3),len(nz)))
        mags[tuple(sorted(np.round(np.abs(ch)**2*e*e).astype(int)))]+=1
    supp[tuple(r[3] for r in row)]+=1
    if shown<maxshow:
        shown+=1; print('SOL',l.split()[1:]); 
        for X,ex in enumerate(E):
            c=[None]*e
            for o in range(1,n): c[cos[o]]=ex[o]
            print('  ',X,'u0=i^%d'%ex[0],'f(j)=',''.join(map(str,c)))
print('v',v,'k',k,'e',e,'g',g)
print('support-size patterns', supp.most_common(8))
print('|e*chat|^2 multisets', mags.most_common(10))

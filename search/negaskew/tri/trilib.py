"""Shared helpers of the tools of this directory that work on the cosets (search/NOTES.md section 92): v = 1 + k(k+1)/2 prime, H = subgroup of order k, e = (v-1)/k cosets.
Coset j = g^j H (g = least primitive root).  A block is (u0, c_0..c_{e-1}) with exponents of i.
load_sols reads the SOL lines of ../qskew (digits per H-orbit, orbits ordered by their smallest element) and of tri2 (digits u0, c_0..c_{e-1} in coset order; the format is recognised by the field 'i=')."""
import sys, numpy as np, itertools, collections, math
from sympy import primitive_root
I4 = np.array([1,1j,-1,-1j])
def setup(v,k):
    g=primitive_root(v); e=(v-1)//k
    assert k*e==v-1
    h=pow(g,e,v); H=[pow(h,i,v) for i in range(k)]
    cos=np.full(v,-1,dtype=int)       # coset index of x != 0 ; cos[0] = -1
    for j in range(e):
        r=pow(g,j,v)
        for x in H: cos[(r*x)%v]=j
    return dict(v=v,k=k,e=e,g=g,H=H,cos=cos)
def load_sols(fn,v,k,orbit_labelling='qdrv'):
    """read qskew SOL lines (orbit exponents, orbits ordered by smallest element, orbit 0 = {0}); return list of 4 blocks each (u0exp, tuple(c exps by coset index j))"""
    d=setup(v,k); e=d['e']; cos=d['cos']
    reps=[]; seen=set()
    orbj=[]
    for x in range(1,v):
        j=cos[x]
        if j not in seen: seen.add(j); orbj.append(j)     # orbits 1..e in order of smallest element
    assert len(orbj)==e
    sols=[]
    for l in open(fn):
        if not l.startswith('SOL'): continue
        P=l.split()[1:]
        if P and P[0].startswith('i='):      # SOL line of tri2: "SOL i=I r=R w=W A B C D", digits u0, c_0..c_{e-1} (coset order)
            sols.append(tuple((int(p[0]),tuple(int(c) for c in p[1:])) for p in P[3:7])); continue
        blk=[]
        for p in P:
            ex=[int(c) for c in p]
            c=[0]*e
            for o in range(1,e+1): c[orbj[o-1]]=ex[o]
            blk.append((ex[0],tuple(c)))
        sols.append(tuple(blk))
    return d,sols
def expand(d,blk):
    v=d['v']; u0,c=blk
    u=np.zeros(v,dtype=complex); u[0]=I4[u0]
    cc=I4[np.array(c)]
    u[1:]=cc[d['cos'][1:]]
    return u
def spectrum(d,blk):
    """Uhat(t) for t = 0, g^0..g^{e-1}"""
    u=expand(d,blk); v=d['v']; g=d['g']; e=d['e']
    F=np.fft.fft(u)           # F[t] = sum u(x) exp(-2 pi i t x / v)  -> equals Uhat with zeta^{-1}; fine (coset relabel)
    return F[0], np.array([F[pow(g,m,v)] for m in range(e)])
def paf(d,blk):
    u=expand(d,blk); v=d['v']
    F=np.fft.fft(u); P=np.fft.ifft(np.abs(F)**2)   # P[s]=sum_x u(x+s) conj u(x)?  (up to sign convention of s)
    return P
def isunit(z): return abs(z)<1e-9

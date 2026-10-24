"""The residue of the cluster at v = 379 under the hash of tri_join.c: reads the first SOL line of a tri2 run (blocks in coset order), computes its four PAF vectors with a slow direct sum over Z_v (the sum must be 4v at 0), and the left residue r = pi(A) + pi(B) mod M, M = 2^18, which equals -(pi(C) + pi(D)); writes runs/r_sol.txt.  This is how the residue 170258 of tri_join 379 27 2 262144 170258 170259 (README.md) was found from a solution of tri2.
usage: python3 test_trijoin379.py [FILE]        (FILE default data/s379.out)"""
import sys, os, numpy as np; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trilib import *
v,k=379,27; d=setup(v,k); e=d['e']; g=d['g']; cos=d['cos']; np_=e//2
# first SOL line of the run
HERE=os.path.dirname(os.path.abspath(__file__)); fn=sys.argv[1] if len(sys.argv)>1 else os.path.join(HERE,'data','s379.out')
for l in open(fn):
    if l.startswith('SOL'): p=l.split(); break
blocks=[]
for s in p[4:8]:
    x=[int(c) for c in s]; blocks.append((x[0],tuple(x[1:])))
def pafvec(b):
    u=expand(d,b); out=[]
    for t in range(np_):
        s=pow(g,t,v); P=sum(u[x]*np.conj(u[(x+s)%v]) for x in range(v)); out+=[int(round(P.real)),int(round(P.imag))]
    return out
Pv=[pafvec(b) for b in blocks]; print('sum of PAF vectors',np.sum(Pv,axis=0))
M64=(1<<64)-1
s=88172645463325252; KR=[];QR=[]
for i in range(128):
    s^=(s<<13)&M64; s^=s>>7; s^=(s<<17)&M64; KR.append(s|1)
    s^=(s<<13)&M64; s^=s>>7; s^=(s<<17)&M64; QR.append(s)
def pi(vec,M): 
    q=sum((x&M64)*QR[i] for i,x in enumerate(vec))&M64; return q&(M-1)
for M in (262144,):
    r=(pi(Pv[0],M)+pi(Pv[1],M))&(M-1); print('M',M,'left residue r =',r,' right residue (pi C + pi D) =',(pi(Pv[2],M)+pi(Pv[3],M))&(M-1),' -r =',(-r)&(M-1))
print('r',r)
os.makedirs(os.path.join(HERE,'runs'),exist_ok=True); open(os.path.join(HERE,'runs','r_sol.txt'),'w').write(str(r))

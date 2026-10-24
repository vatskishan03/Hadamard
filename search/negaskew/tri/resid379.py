"""Distribution of the residues of the hit partitions: for the first 300 solutions of a tri2 run at v = 379 the left residue r = pi(A) + pi(B) mod 2^18 under the hash of tri_join.c (xorshift seed 88172645463325252), and its distribution mod 2, 4, ..., 64 (search/NOTES.md section 92: the solutions are spread over the residues almost at random).
usage: python3 resid379.py [FILE]        (FILE default data/s379.out)"""
import sys, os, numpy as np, collections; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trilib import *
v,k=379,27; d=setup(v,k); e=d['e']; g=d['g']; np_=e//2
M64=(1<<64)-1
s=88172645463325252; KR=[];QR=[]
for i in range(128):
    s^=(s<<13)&M64; s^=s>>7; s^=(s<<17)&M64; KR.append(s|1)
    s^=(s<<13)&M64; s^=s>>7; s^=(s<<17)&M64; QR.append(s)
# precompute N tensors for fast PAF
orb=np.zeros(v,dtype=np.int64)
for x in range(1,v): orb[x]=d['cos'][x]+1
Nt=np.zeros((np_,e+1,e+1))
for t in range(np_):
    s_=pow(g,t,v); np.add.at(Nt[t],(orb,orb[(np.arange(v)+s_)%v]),1)
def pafv(b):
    ex=np.zeros(e+1,dtype=np.int64); ex[0]=b[0]; ex[1:]=b[1]; U=I4[ex]
    out=[]
    for t in range(np_):
        P=(U@Nt[t]*np.conj(U)).sum()
        # the convention of tri_join.c: sum u(x) conj(u(x+s)) with a=orbit(x), b=orbit(x+s): N[a][b] u_a conj(u_b)
        out+=[int(round(P.real)),int(round(P.imag))]
    return out
def pi(vec,M): 
    q=sum((x&M64)*QR[i] for i,x in enumerate(vec))&M64; return q&(M-1)
M=262144; rs=[]
fn=sys.argv[1] if len(sys.argv)>1 else os.path.join(os.path.dirname(os.path.abspath(__file__)),'data','s379.out')
for l in open(fn):
    if not l.startswith('SOL'): continue
    p=l.split(); bl=[]
    for s_ in p[4:8]:
        x=[int(c) for c in s_]; bl.append((x[0],tuple(x[1:])))
    Pv=[pafv(b) for b in bl]
    rs.append((pi(Pv[0],M)+pi(Pv[1],M))&(M-1))
    if len(rs)>=300: break
rs=np.array(rs); print('n',len(rs))
for kb in (1,2,3,4,5,6):
    c=collections.Counter(rs%(1<<kb)); print('r mod',1<<kb,sorted(c.items()))
print('min r',rs.min(),'quantiles',np.quantile(rs,[0.1,0.5,0.9]))

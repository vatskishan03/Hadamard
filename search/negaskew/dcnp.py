"""Row-sum (DC) feasibility of an orbit shape for the quaternary skew problem of README.md.  A block constant on the
H-orbits has row sum R = sum_a |O_a| u_a with u_a a unit i^e of Z[i] (the skew block A: R_A = u(0) + (1+i) * (integers),
u(0) in {1, -i}); the equation sum_X |R_X|^2 = 4v needs a combination of one R_A and three R_B/R_C/R_D values.
plain_dc(v, sz, neg, n)  True iff some (R_A, R_B, R_C, R_D) has sum |R|^2 = 4v (an EXACT necessary condition);
spin_dc(v, sz, neg, n)   the same for the spin structure C = B o mu, D = B o mu^2 (|R_A|^2 + 3 |R_B|^2 = 4v).
usage: python3 dcnp.py     prints plain_dc for v = 89 <78> (False: no row-sum pattern), 191 <180> and 277 <273> (True)."""
import numpy as np
def rset_B(v,sz,n):
    S=2*v+1; off=v
    cur=np.zeros((S,S),dtype=bool); cur[off+sz[0],off]=True   # w0=1 -> (sz0,0)
    for a in range(1,n):
        s=sz[a]; nxt=np.zeros_like(cur)
        for dx,dy in ((s,0),(0,s),(-s,0),(0,-s)):
            nxt|=np.roll(np.roll(cur,dx,axis=0),dy,axis=1)
        cur=nxt
    xs,ys=np.nonzero(cur); return set(zip((xs-off).tolist(),(ys-off).tolist()))
def rset_A(v,sz,neg,n):
    S=2*v+1; off=v
    cur=np.zeros((S,S),dtype=bool); cur[off+1,off]=True; cur[off,off-1]=True
    for a in range(1,n):
        if a<neg[a]:
            s=sz[a]; nxt=np.zeros_like(cur)
            for d in (s,-s):
                nxt|=np.roll(np.roll(cur,d,axis=0),d,axis=1)
            cur=nxt
    xs,ys=np.nonzero(cur); return set(zip((xs-off).tolist(),(ys-off).tolist()))
def spin_dc(v,sz,neg,n):
    RB=rset_B(v,sz,n); RA=rset_A(v,sz,neg,n)
    nA={x*x+y*y for x,y in RA}
    tot=4*v
    for x,y in RB:
        r=tot-3*(x*x+y*y)
        if r>=0 and r in nA: return True
    return False
def plain_dc(v,sz,neg,n):
    RB=rset_B(v,sz,n); RA=rset_A(v,sz,neg,n)
    nB={}
    for x,y in RB: nB[x*x+y*y]=1
    nA={x*x+y*y for x,y in RA}
    tot=4*v
    ks=sorted(nB)
    two=set()
    for a in ks:
        if a>tot: break
        for b in ks:
            if a+b>tot: break
            two.add(a+b)
    for na in nA:
        r=tot-na
        if r<0: continue
        for t in two:
            if r-t in nB: return True
    return False
if __name__=="__main__":
    from qdrv import shape
    for v,g in [(89,[78]),(191,[180]),(277,[273])]:
        d=shape(v,g); print(v,plain_dc(v,d['sz'],d['neg'],d['n']))

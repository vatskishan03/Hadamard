"""Exact cost of the plain four-list join of qskew for an orbit shape: the number of table entries (build) and of probes,
from the distributions of the row sums (a dynamic programme over the orbit sizes), with the row-sum buckets of qskew.
usage: python3 costest.py v:[gens] ...      e.g.  python3 costest.py 191:[180] 277:[273] 299:[289]
Prints n, build, probe (and their log2), |L| = size of a free list, |A| = size of the skew list.  The probes are about
8e9 at 191 (n = 11), 7.9e12 at 277 (n = 13, where the hit came after about 1% of the work) and 6.8e15 at 299 (n = 15):
n = 15 is some 10^3 times n = 13, out of reach of the join (NOTES.md section 88)."""
import numpy as np, sys, json
from qdrv import shape
def dpB(v,sz,n):
    S=2*v+1; off=v
    cur=np.zeros((S,S),dtype=np.float64); cur[off+sz[0],off]=1.0
    for a in range(1,n):
        s=sz[a]; nxt=np.zeros_like(cur)
        for dx,dy in ((s,0),(0,s),(-s,0),(0,-s)):
            nxt+=np.roll(np.roll(cur,dx,axis=0),dy,axis=1)
        cur=nxt
    return cur
def dpA(v,sz,neg,n):
    S=2*v+1; off=v
    cur=np.zeros((S,S),dtype=np.float64); cur[off+1,off]=1; cur[off,off-1]=1
    for a in range(1,n):
        if a<neg[a]:
            s=sz[a]; nxt=np.zeros_like(cur)
            # e_a=j -> unit i^j ; partner e_neg = 1-j -> i^(1-j): sum = s*(i^j + i^(1-j))
            for j in range(4):
                z=(1j)**j+(1j)**((1-j)%4); dx,dy=int(round(z.real))*s,int(round(z.imag))*s
                nxt+=np.roll(np.roll(cur,dx,axis=0),dy,axis=1)
            cur=nxt
    return cur
def norm_hist(cnt,v):
    S=2*v+1; off=v
    xs,ys=np.nonzero(cnt)
    h={}
    for x,y in zip(xs,ys):
        nm=(x-off)**2+(y-off)**2
        h[nm]=h.get(nm,0)+cnt[x,y]
    return h
def cost(v,g):
    d=shape(v,g); n=d['n']; sz=d['sz']; neg=d['neg']
    cB=dpB(v,sz,n); cA=dpA(v,sz,neg,n)
    hB=norm_hist(cB,v); hA=norm_hist(cA,v)
    tot=4*v
    # counts of R-values per norm (number of buckets) not needed; probe cost uses element counts
    # probes for pairing (A B | C D): sum over (RA,RB) buckets,(RC,RD): need per-bucket sizes; approximate by norm-aggregated sums
    # number of buckets with given norm needed for table multiplicity: build side sum_{valid (RA,RB)} |A||B| = sum over norms hA[na]*hB[nb] (aggregated over all buckets with these norms)
    # probe side: for each table bucket (RA,RB) [norms na,nb], probes = sum_{nc,nd: sum ok} hB[nc]*hB[nd] -> total = sum_{na,nb} (#bucket pairs with norms na,nb) * ...
    # get number of buckets per norm
    def nb(cnt):
        xs,ys=np.nonzero(cnt); h={}
        for x,y in zip(xs,ys):
            nm=(x-v)**2+(y-v)**2; h[nm]=h.get(nm,0)+1
        return h
    bB=nb(cB); bA=nb(cA)
    # T2[s] = sum_{nc+nd=s} hB[nc]hB[nd] (element pairs) ; T2b[s] = sum bB bB (bucket pairs)
    def conv(h1,h2,lim):
        out={}
        for a,x in h1.items():
            if a>lim: continue
            for b,y in h2.items():
                if a+b>lim: continue
                out[a+b]=out.get(a+b,0)+x*y
        return out
    P_AB_elem=conv(hA,hB,tot); P_AB_bk=conv(bA,bB,tot)
    P_CD_elem=conv(hB,hB,tot)
    build=0; probe=0
    for s,eb in P_AB_elem.items():
        if tot-s in P_CD_elem:
            build+=eb
            # probes: for each (RA,RB) bucket pair with norm sum s: sum of |C||D| over (RC,RD) with norm sum tot-s => (#bucket pairs)*P_CD_elem[tot-s]
            probe+=P_AB_bk[s]*P_CD_elem[tot-s]
    return n,build,probe,sum(hB.values()),sum(hA.values())
if __name__=="__main__":
    for arg in sys.argv[1:]:
        v,g=arg.split(':'); v=int(v); g=json.loads(g)
        n,b,p,nl,na=cost(v,g)
        print(v,g,'n',n,'build %.3g probe %.3g (log2 %.1f %.1f) |L|=%.3g |A|=%.3g'%(b,p,np.log2(max(b,1)),np.log2(max(p,1)),nl,na),flush=True)

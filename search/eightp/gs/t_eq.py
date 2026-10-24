import numpy as np
from gsk import *
rng=np.random.default_rng(1)
for (Km,r) in [([5],3),([4],5),([3,3],7),([6],9),([2,2],3),([7],13)]:
    sh=GShape(Km,r)
    N=sh.N
    bad=0
    for trial in range(20):
        blocks=[]
        for b in range(4):
            x0=rng.choice([-1,1],N); x1=rng.choice([-1,1],N); x2=rng.choice([-1,1],N)
            e=(x1+x2)//2; t=(x1-x2)//2; s=x0-e
            blocks.append((s,t,e))
        seqs=sh.glevel_blocks(blocks)
        paf=exact_paf_G(sh.gmods(),seqs)
        U,I,J,Q=[sum(sh.full_eqs(*b)[i] for b in blocks) for i in range(4)]
        p,ee=factor_prime_power(r)
        F=GF(p,ee); sq={int(F.mul[x,x]) for x in range(1,r)}
        pred=np.zeros(N*r,dtype=np.int64)
        for k in range(N):
            for a in range(r):
                chi=0 if a==0 else (1 if a in sq else -1)
                pred[k*r+a]=I[k]*(a==0)+Q[k]*chi+J[k]
        # alt convention Q(-k)
        pred2=np.zeros(N*r,dtype=np.int64)
        for k in range(N):
            for a in range(r):
                chi=0 if a==0 else (1 if a in sq else -1)
                pred2[k*r+a]=I[k]*(a==0)+Q[sh.negperm[k]]*chi+J[k]
        if not (pred==paf).all(): bad+=1
        if not (U==I+r*J).all(): bad+=100
    print(Km,r,'mismatches',bad)

# verify the defect identity (r+1) e_p + (r-3) e_m = 4|K| + phi(0) on actual solutions (K-level, from the exhaustive engine) and the row-sum identities
import sys; sys.path.insert(0,'.')
from gsk import *; from gs_exh import exhaust
import numpy as np
cases=[('spence',[4],3,[],None),('spence',[2,2],3,[],None),
       ('bord',[3],3,[],[(0,),(1,),(2,)]),('bord',[2,2],3,[],[(1,0),(0,1),(1,1)]),('bord',[6],5,[],[(0,),(2,),(4,)]),('bord',[4],3,[],[(1,),(3,)])]
for (mode,K,r,hg,T) in cases:
    sh=GShape(K,r,hg,mode,T)
    tot,sols=exhaust(sh,want_sols=300,tmax=40,max_inst=100000)
    bad=0; n=0
    for bl in sols:
        ep=em=0; rs=[]
        for (s,t,e) in bl:
            d=(e!=0); ep+=int((d&(s==0)).sum()); em+=int((d&(s!=0)).sum())
        if (r+1)*ep+(r-3)*em!=4*sh.N+int(sh.phi[0]): bad+=1
        # row sums of the four G-level sequences: sum R^2 == Rsq
        seqs=sh.glevel_blocks(bl); R=[int(x.sum()) for x in seqs]
        if sum(x*x for x in R)!=sh.Rsq: bad+=1000
        n+=1
    print(mode,K,r,hg,T,'phi0',int(sh.phi[0]),'solutions checked',n,'identity/row-sum violations',bad)

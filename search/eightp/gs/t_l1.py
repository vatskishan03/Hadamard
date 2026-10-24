# test level-1: explicit lists all four blocks vs numpy exact pair-join count on tiny shapes
import numpy as np, itertools, sys, time
from gsk import *
from gl1 import *
def run(Km,r,hg=(),mode='gs',T=None,maxal=60):
    sh=GShape(Km,r,hg,mode,T); l1=L1(sh)
    lists={}
    for n in range(0,sh.N+1):
        rr=enum_weight(sh,n)
        if rr is not None: lists[n]=rr
    def feats(S,E):
        u=S.astype(np.int64)+sh.r*E.astype(np.int64)
        return np.stack([(u*u[:,sh.shift[k]]).sum(1) for k in sh.repsPM],1)
    def enc(F): 
        off=1<<14
        out=np.zeros(len(F),np.int64)
        for j in range(F.shape[1]): out=out*(1<<16)+(F[:,j]+off)
        return out
    tot=0; nal=0; allok=True
    for al in itertools.combinations_with_replacement(sorted(lists.keys()),4):
        if sum(al) not in sh.feasible_totals(): continue
        specs=[BlockSpec(al[b],S=lists[al[b]][0],E=lists[al[b]][1]) for b in range(4)]
        cnt,S,E,st=l1.run(specs,maxout=1)
        if max(len(lists[al[0]][0])*len(lists[al[2]][0]),len(lists[al[1]][0])*len(lists[al[3]][0]))>2e6: continue
        F=[feats(*lists[a]) for a in al]
        # exact: encode pair sums; need sum==Ut: enc is not additive, so use unique rows
        L=(F[0][:,None,:]+F[2][None,:,:]).reshape(-1,sh.nU); Rr=(F[1][:,None,:]+F[3][None,:,:]).reshape(-1,sh.nU)
        need=sh.Utarget[None,:]-L
        eL=enc(need); eR=enc(Rr)
        u,c=np.unique(eR,return_counts=True)
        pos=np.searchsorted(u,eL); pos[pos>=len(u)]=len(u)-1
        ex=int(c[pos][u[pos]==eL].sum())
        ok=(cnt==ex); allok&=ok; tot+=cnt; nal+=1
        if nal<=maxal and (cnt or not ok): print(' ',al,cnt,ex,'OK' if ok else 'MISMATCH')
    print(Km,r,hg,mode,'allocs',nal,'total level-1 instances',tot,'ALLOK' if allok else 'MISMATCH')
import os
if os.environ.get('FULL') == '1':
    run([4],3); run([3],5); run([5],3,maxal=5)
else:
    run([3],5); run([2,2],3,mode='spence'); run([3],3,mode='bord',T=[(0,),(1,),(2,)]); run([2],5)

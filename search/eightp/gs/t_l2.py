# level-2 test: level-1 instances -> gl2 completions vs vectorised brute force over all t-sign assignments
import numpy as np, itertools, sys, time
from gsk import *
from gl1 import *
from gl2io import *
rng=np.random.default_rng(7)
def bf_count(sh,blocks,maxbits=15):
    K=sh.kCD
    free=[[i for i in range(K.nu) if not (e[int(K.pos[i][0])]!=0)] for (s,e) in blocks]
    nv=[len(f) for f in free]; tot=sum(nv)
    if tot>maxbits: return None
    M=1<<tot; ids=np.arange(M)
    N=sh.N; r=sh.r; eps=sh.eps
    Isum=np.zeros((M,N),np.int64); Qsum=np.zeros((M,N),np.int64); Jsum=np.zeros((M,N),np.int64)
    pos=0
    for b,(s,e) in enumerate(blocks):
        t=np.zeros((M,N),np.int64)
        for j,i in enumerate(free[b]):
            P=K.pos[i]; bit=1-2*((ids>>(pos+j))&1)
            t[:,P]=bit[:,None]*K.sig[P][None,:]
        pos+=nv[b]
        s64=s.astype(np.int64); e64=e.astype(np.int64)
        for k in range(N):
            sh_=sh.shift[k]
            paf=lambda x: (x*x[:,sh_]).sum(1)
            Isum[:,k]+= (s64*s64[sh_]).sum()+r*paf(t)
            Jsum[:,k]+= r*(e64*e64[sh_]).sum()+(s64*e64[sh_]).sum()+(e64*s64[sh_]).sum()-paf(t)
            Qsum[:,k]+= (s64[None,:]*t[:,sh_]).sum(1)+eps*(t*s64[sh_][None,:]).sum(1)   # ccf(s,t)+eps*ccf(t,s): sum_j s(j)t(j+k)+eps sum_j t(j)s(j+k)
    tI=np.zeros(N,np.int64); tI[0]=sh.T2
    ok=(Isum==tI[None,:]).all(1)&(Qsum==0).all(1)&(Jsum==sh.phi[None,:]).all(1)
    return int(ok.sum())
def test(Km,r,hg,mode,T,ninst=40,maxbits=15):
    sh=GShape(Km,r,hg,mode,T); l1=L1(sh)
    lists={n:enum_weight(sh,n) for n in range(sh.N+1)}
    lists={n:v for n,v in lists.items() if v is not None}
    insts=[]
    for al in itertools.combinations_with_replacement(sorted(lists),4):
        if sum(al) not in sh.feasible_totals(): continue
        specs=[BlockSpec(al[b],S=lists[al[b]][0],E=lists[al[b]][1]) for b in range(4)]
        if max(specs[0].L*specs[2].L, specs[1].L*specs[3].L)>3e5: continue
        cnt,S,E,st=l1.run(specs,maxout=3000)
        for i in range(min(cnt,len(S))):
            if rng.random()<0.1: insts.append([(S[i,b],E[i,b]) for b in range(4)])
        if len(insts)>=ninst*3: break
    # prefer instances with few free variables so that brute force is possible
    K=sh.kCD
    def nvars(bl): return sum(int(sum(1 for i in range(K.nu) if e[int(K.pos[i][0])]==0)) for (s,e) in bl)
    insts=[b for b in insts if nvars(b)<=maxbits][:ninst]
    res,sols=run_gl2(sh,insts,maxsol=100)
    nm=0; tot=0; nsol=0; nbf=0
    for i,bl in enumerate(insts):
        c=res[i]['compl']; tot+=c
        b=bf_count(sh,[(s,e) for (s,e) in bl],maxbits)
        nbf+=1
        # gl2 reports completions as t-quadruples counted with the global t->-t gauge ALREADY multiplied (x2)
        if b!=c: nm+=1; print('MISMATCH inst',i,'bf',b,'gl2',c,res[i]['raw'])
    for (i,ts) in sols:
        bl=[(insts[i][X][0].astype(np.int64),ts[X].astype(np.int64),insts[i][X][1].astype(np.int64)) for X in range(4)]
        assert sh.check_klevel(bl), 'solution fails K-level check'
        nsol+=1
    print(Km,r,hg,mode,'insts',len(insts),'compared',nbf,'mismatch',nm,'completions',tot,'sols verified',nsol,flush=True)
if __name__=='__main__':
    t0=time.time()
    test([5],3,[],'gs',None)
    test([4],3,[],'gs',None)
    test([2,2],3,[],'spence',None)
    test([3],5,[],'gs',None)
    test([7],3,[([2],1)],'gs',None)
    print(time.time()-t0)

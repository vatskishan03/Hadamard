"""Defect patterns of the Paley plug-in in the quaternary skew problem on Z_v = Z_K x F_q, H = {1} x QR(q) x H_K: all
defect configurations (K-orbit, type) of the blocks with total weight 4|K| (the defect identity), up to dmax defects per block.
usage: python3 alloc.py v K q "[gens]" dmax       (gens: generators of H in Z_v, e.g. from qshape.cyc_gens)
Prints the number of K-orbits, the orbit sizes and the number of defect-configuration quadruples of weight sum 4|K|;
0 means no solution at all (EXACTLY empty: 295, 335, 515, 655, 695, 995, 927, 611, 899 of NOTES.md section 89).
list_patterns(v, K, q, gens, dmax) is used by pal_lists.py."""
import sys, json, itertools
from math import comb
from idfeas import korbits
from enum_shapes import gen_sub
def allocs(v,K,q,HK,maxprint=40):
    orbs,idx=korbits(K,HK); kap=len(orbs); sz=[len(o) for o in orbs]
    neg=[idx[(-o[0])%K] for o in orbs]
    total=4*K
    Wf=[('F+',q+1),('F-',q-3),('Fi',q-1)]; Wm=[('M+',(q+1)//2),('M-',(q-3)//2)]
    # block signature = tuple of (orbit index, type) ; enumerate per-block defect sets with weights, group by (multiset of (size,weight))
    # per block B/C/D: choose a subset of orbits and a type for each
    opts=[(n,w) for n,w in Wf+Wm]
    # generate all per-block defect configurations up to dmax defects with weight <= total
    def block_cfgs(dmax,paired):
        res=[]  # list of (weight, tuple of (orbit,type))
        orbs_list=[i for i in range(kap) if (not paired) or (i>0 and i<neg[i])]
        def rec(pos,cur,wt):
            res.append((wt,tuple(cur)))
            if len(cur)>=dmax: return
            for j in range(pos,len(orbs_list)):
                i=orbs_list[j]
                for nme,w in opts:
                    c=(2 if paired else 1)*sz[i]*w
                    if wt+c<=total:
                        cur.append((i,nme)); rec(j+1,cur,wt+c); cur.pop()
        rec(0,[],0)
        return res
    return orbs,sz,block_cfgs
if __name__=="__main__":
    v=int(sys.argv[1]); K=int(sys.argv[2]); q=int(sys.argv[3]); g=json.loads(sys.argv[4]); dmax=int(sys.argv[5])
    H=gen_sub(v,g); HK=sorted(set(h%K for h in H))
    orbs,sz,bc=allocs(v,K,q,HK)
    print('kappa',len(orbs),'sizes',sz)
    cfgB=bc(dmax,False); cfgA=bc(dmax,True)
    # A: forced @0 mixed weights
    Wm=[(q+1)//2,(q-3)//2]
    from collections import defaultdict
    # weight histogram
    hB=defaultdict(int)
    for w,c in cfgB: hB[w]+=1
    hA=defaultdict(int)
    for w,c in cfgA:
        for w0 in Wm: hA[w+w0]+=1
    total=4*K
    # count weight-compatible quadruples (as multiplicities of configurations)
    # convolve B,C,D
    import numpy as np
    tot=total
    def hist_arr(h):
        a=np.zeros(tot+1,dtype=object)
        for w,c in h.items():
            if w<=tot: a[w]+=c
        return a
    a=hist_arr(hA); b=hist_arr(hB)
    def conv(x,y):
        z=np.zeros(tot+1,dtype=object)
        for i in range(tot+1):
            if x[i]:
                for j in range(tot+1-i):
                    if y[j]: z[i+j]+=x[i]*y[j]
        return z
    bb=conv(b,b); bbb=conv(bb,b); allc=conv(a,bbb)
    print('defect-config quadruples with weight sum = 4K:',allc[tot], '  (blocks with <=%d defects)'%dmax)
    print('A cfgs',len(cfgA),'B cfgs',len(cfgB))

def list_patterns(v,K,q,g,dmax):
    H=gen_sub(v,g); HK=sorted(set(h%K for h in H))
    orbs,sz,bc=allocs(v,K,q,HK)
    cfgB=bc(dmax,False); cfgA=bc(dmax,True)
    Wm=[(q+1)//2,(q-3)//2]; tot=4*K
    from collections import defaultdict
    byw=defaultdict(list)
    for w,c in cfgB: byw[w].append(c)
    pats=[]
    for wA,cA in cfgA:
        for w0 in Wm:
            for wB,cB_ in cfgB:
                for wC,cC in cfgB:
                    wD=tot-wA-w0-wB-wC
                    if wD<0 or wD not in byw: continue
                    for cD in byw[wD]:
                        pats.append(((cA,w0),cB_,cC,cD))
    return pats,sz

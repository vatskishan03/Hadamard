"""check the unit-action permutations at the tail of shape files written by qshape37.py / qshape37m.py: every permutation must be a bijection of the K-orbits that preserves the orbit sizes, and the set (with the identity) must be closed under composition.
usage: python3 unitcheck.py 'shapes/*.txt'   (the check that found the uA defect of the first qshape37m.py: all diagonal unit pairs, not only those normalising H_K)"""
import sys, glob
def parse(fn):
    toks=open(fn).read().split()
    nu=int(toks[0]); q=int(toks[1]); K=int(toks[2]); pos=3
    ell=list(map(int,toks[pos:pos+nu])); pos+=nu
    neg=list(map(int,toks[pos:pos+nu])); pos+=nu
    pos+=nu**3  # N tensor
    pos+=2*nu*nu  # C matrix (re,im)
    if pos>=len(toks): return nu,q,K,ell,neg,None
    nU=int(toks[pos]); pos+=1
    U=[tuple(map(int,toks[pos+i*nu:pos+(i+1)*nu])) for i in range(nU)]
    return nu,q,K,ell,neg,U
bad=[]; tot=0
for fn in sorted(glob.glob(sys.argv[1])):
    try: nu,q,K,ell,neg,U=parse(fn)
    except Exception as e: continue
    if not U: continue
    tot+=1
    ident=tuple(range(nu)); S=set(U)|{ident}
    prob=[]
    for p in U:
        if sorted(p)!=list(range(nu)): prob.append('notbij')
        elif [ell[p[o]] for o in range(nu)]!=ell: prob.append('ell')
    if not prob:
        # closure
        for a in U:
            for b in U:
                c=tuple(a[b[i]] for i in range(nu))
                if c not in S: prob.append('notclosed'); break
            if prob: break
    if prob: bad.append((fn,len(U),sorted(set(prob))))
print('checked',tot,'shapes with unit perms; invalid:',len(bad))
for b in bad: print(b)

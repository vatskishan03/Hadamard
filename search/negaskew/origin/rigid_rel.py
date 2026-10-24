"""relations among blocks in rigid (A,M|E,E)-type solutions: python3 rigid_rel.py m q "[hK]" solfile
 tests, for each ordered pair of 'functions' among s_X, t_X (off the origin) the relation f = rho * g^{(*)}(sigma k)  with sigma in {+k, -k}, (*) in {none, conj}, rho unit constant"""
import sys, itertools
sys.path.insert(0,'.')
from common import *
from collections import Counter
m=int(sys.argv[1]); q=int(sys.argv[2]); hK=eval(sys.argv[3]); fn=sys.argv[4]
S=Shape(m,hK); orb=np.array(S.orb)
sols=load_solutions(fn)
import os
sols=sols[::int(os.environ.get("STEP","1"))]
print(len(sols),'solutions')
Uarr=np.array(U)
def typ(b):
    d=(b[2][0]-b[1][0])%4
    e=(U[b[1][0]]+U[b[2][0]])/2; x=U[b[0][0]]
    c=round((q-1)*abs(e)**2+2*(x*e.conjugate()).real)
    return {2:'E',0:'F',1:'M',3:'M'}[d]+str(c)
pat=Counter(); rel=Counter()
neg=[(-k)%m for k in range(m)]
for sol in sols:
    pat[tuple(typ(b) for b in sol)]+=1
    fs={}
    for X,b in enumerate(sol):
        x0=Uarr[np.array(b[0])[orb]]; x1=Uarr[np.array(b[1])[orb]]; x2=Uarr[np.array(b[2])[orb]]
        e=(x1+x2)/2; s=x0-e; t=(x1-x2)/2
        fs['s'+'ABCD'[X]]=s[1:]; fs['t'+'ABCD'[X]]=t[1:]
    names=list(fs)
    for f,g in itertools.combinations(names,2):
        for refl in (0,1):
            for cj in (0,1):
                gg=fs[g][np.array(neg[1:])-1] if refl else fs[g]
                gg=np.conj(gg) if cj else gg
                ratio=fs[f]/gg
                if np.max(np.abs(ratio-ratio[0]))<1e-9 and abs(abs(ratio[0])-1)<1e-9:
                    rel[(f,g,refl,cj)]+=1
print(pat.most_common(5))
tot=len(sols)
for k,v in rel.most_common(25): print(k,v,round(v/tot,3))

import sys, itertools
sys.path.insert(0,'.')
from common import *
from collections import Counter
m=19; S=Shape(19,[7]); orb=S.orb; nu=S.nu; neg=S.sh['neg']
from sympy import primitive_root
gr=primitive_root(m); e=nu-1
def ind(k):
    x=1
    for j in range(m-1):
        if x==k: return j
        x=x*gr%m
cidx={orb[k]:ind(k)%e for k in range(1,m)}
inv={c:o for o,c in cidx.items()}
sols=load_solutions('runs/589_sample.out'); S0=set(sols)
def mult(sol,s):
    out=[]
    for b in sol:
        nb=[]
        for comp in b:
            new=list(comp)
            for o in range(1,nu): new[o]=comp[inv[(cidx[o]+s)%e]]
            nb.append(tuple(new))
        out.append(tuple(nb))
    return tuple(out)
def swap12(sol): return tuple((b[0],b[2],b[1]) for b in sol)
def normalize(sol):
    out=[]
    for X,b in enumerate(sol):
        if X==0: out.append(b); continue
        c=b[0][0]; out.append(tuple(tuple((x-c)%4 for x in comp) for comp in b))
    return tuple(out)
def conjall(sol):
    out=[]
    for X,b in enumerate(sol):
        nb=tuple(tuple((-x)%4 for x in comp) for comp in b)
        if X==0: nb=tuple(tuple((x+3)%4 for x in comp) for comp in nb)
        out.append(nb)
    return normalize(tuple(out))
def reflect(sol):
    return normalize(tuple(tuple(tuple(comp[neg[o]] for o in range(nu)) for comp in b) for b in sol))
def perm(sol,p): return (sol[0],)+tuple(sol[1+i] for i in p)
perms=list(itertools.permutations(range(3)))
def images(s):
    out=set()
    for sh in range(e):
        a=mult(s,sh)
        for sw in range(2):
            b=swap12(a) if sw else a
            for cj in range(2):
                c=conjall(b) if cj else b
                for rf in range(2):
                    d=reflect(c) if rf else c
                    for p in perms: out.add(perm(d,p))
    return out
# group-closure check on the sample members
print('sample',len(sols))
seen=set(); classes=[]
cnt=Counter()
for s in sols:
    if s in seen: continue
    im=images(s)
    inside=[x for x in im if x in S0]
    seen|=set(inside)
    classes.append(len(inside))
    cnt[len(im)]+=1
print('distinct classes in the sample:',len(classes),' class sizes (group orbit sizes):',cnt, ' members of each class present in the sample:',Counter(classes))

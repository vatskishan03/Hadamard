import sys, itertools
sys.path.insert(0,'.')
from common import *
from collections import Counter, defaultdict
m=13;q=23
S=Shape(13,[3])
orb=S.orb; nu=S.nu
neg=S.sh['neg']
sols=load_solutions('runs/299_all.out')
# class index c(k) = log_2 k mod 4 ; orbit o = c+1 ; check
def ind2(k):
    x=1
    for j in range(12):
        if x==k: return j
        x=x*2%13
cidx={}
for k in range(1,13): cidx[orb[k]]=ind2(k)%4
print('orbit->class', cidx)
# multiplier by 2^s : maps orbit of class c to class c+s ; the block X'(k)=X(2k)?? -> X'(orbit c)=X(orbit c+1)
def mult(sol, s):
    # X'(k) = X(mu k), mu=2^s: orbit o (class c) takes data of class c+s
    inv={c:o for o,c in cidx.items()}
    out=[]
    for b in sol:
        nb=[]
        for comp in b:
            new=list(comp)
            for o in range(1,nu):
                new[o]=comp[inv[(cidx[o]+s)%4]]
            nb.append(tuple(new))
        out.append(tuple(nb))
    return tuple(out)
def swap12(sol):
    return tuple((b[0],b[2],b[1]) for b in sol)
def normalize(sol):
    out=[]
    for X,b in enumerate(sol):
        if X==0: out.append(b); continue
        c=b[0][0]
        out.append(tuple(tuple((x-c)%4 for x in comp) for comp in b))
    return tuple(out)
def conjall(sol):
    out=[]
    for X,b in enumerate(sol):
        nb=tuple(tuple((-x)%4 for x in comp) for comp in b)
        if X==0: nb=tuple(tuple((x+3)%4 for x in comp) for comp in nb)
        out.append(nb)
    return normalize(tuple(out))
def reflect(sol):
    out=[]
    for X,b in enumerate(sol):
        nb=tuple(tuple(comp[neg[o]] for o in range(nu)) for comp in b)
        out.append(nb)
    return normalize(tuple(out))
def perm(sol,p):
    return (sol[0],)+tuple(sol[1+i] for i in p)
S0=set(sols)
print(len(S0))
# check closure of generators
def closed(f):
    bad=0
    for s in sols:
        if f(s) not in S0: bad+=1
    return bad
print('mult1 bad',closed(lambda s:mult(s,1)))
print('swap bad',closed(swap12))
print('conj bad',closed(conjall))
print('reflect bad',closed(reflect))
print('perm bad',closed(lambda s:perm(s,(1,0,2))), closed(lambda s:perm(s,(0,2,1))))

perms=list(itertools.permutations(range(3)))
def images(s):
    out=set()
    for sh in range(4):
        a=mult(s,sh)
        for sw in range(2):
            b=swap12(a) if sw else a
            for cj in range(2):
                c=conjall(b) if cj else b
                for rf in range(2):
                    d=reflect(c) if rf else c
                    for p in perms:
                        out.add(perm(d,p))
    return out
seen=set(); classes=[]
for s in sols:
    if s in seen: continue
    im=images(s)
    seen|=im
    classes.append((s,len(im)))
print('classes',len(classes), Counter(n for _,n in classes))
import pickle
pickle.dump(classes,open('runs/classes299.pkl','wb'))
for s,n in classes[:5]: print(n,s)

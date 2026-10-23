"""Classification of the multiplier-type symmetries a(x) = eps x^s a(x^t) of a palindromic (skew) negacyclic sequence a of length 4v (a_0 = 1): count_inv(v, gens)
is the number of palindromic a invariant under all (eps, s, t) of gens.  usage: python3 mult_classify.py     for v = 3, 5, 7, 9, 15, 21, 35 the elements (eps, s, t) other than
the identity that admit an invariant skew sequence, by (t mod 8, t mod v): only t = 1 (mod 8) admits many (>= 16); the Galois part t = 3, 5, 7 (mod 8) never does
(NOTES.md section 89: H <= Z_v* embedded with t = 1 (mod 8), unshifted, -1 not in H)."""
import sys, itertools
from math import gcd
sys.setrecursionlimit(10000)
def count_inv(v, gens):
    """gens: list of (eps,s,t). Return number of palindromic a (a0=1) with a = eps x^s a(x^t) for all gens (or 0).
    a viewed as antiperiodic function on Z_{8v}: a[j+4v] = -a[j]. Palindromic: a_j = a_{4v-j} (j=1..4v-1), a_0=1."""
    N=8*v; h=4*v
    # unknown a_j for j in 0..N-1 with relation a[j+h] = - a[j]; use union-find with parity on indices j mod N
    parent=list(range(N)); par=[0]*N   # par = sign parity relative to parent (1 = opposite)
    bad=[False]
    def find(x):
        if parent[x]==x: return x,0
        r,p=find(parent[x]); parent[x]=r; par[x]^=p; return r,par[x]
    def union(x,y,sgn): # a[x] = (-1)^sgn a[y]
        rx,px=find(x); ry,py=find(y)
        if rx==ry:
            if (px^py)!=sgn: bad[0]=True
        else:
            parent[rx]=ry; par[rx]=px^py^sgn
    for j in range(N): union((j+h)%N,j,1)
    # palindromic: a_j = a_{4v-j}  for j=1..4v-1 ;  as antiperiodic: a_{-j} = -a_{4v-j} => a_{4v-j} = -a_{-j}: a_j = - a_{-j}
    for j in range(1,h): union(j,(-j)%N,1)
    for (eps,s,t) in gens:
        # a(x) = eps x^s a(x^t): coefficient at x^m: a_m = eps * a_{(m-s)/t}
        ti=pow(t,-1,N)
        for m in range(N):
            j=((m-s)*ti)%N
            union(m,j,0 if eps==1 else 1)
    # a_0 = 1 fixed: need a_0 not forced opposite to itself/ a_{4v}: fine. a_0 and a_{4v} = -a_0 handled. 
    # check a_0 != its own contradictions: bad flag. Also a_0 must be +1: no constraint except sign choice. Count classes
    if bad[0]: return 0
    roots=set(find(j)[0] for j in range(N))
    # classes include both j and j+h in same class (parity): number of free classes = len(roots) (each class free sign), minus 1 for a_0=+1 fixed
    return 2**(len(roots)-1) if len(roots)>=1 else 0
def order_of(v,eps,s,t):
    # order of the map on sequences
    N=8*v
    # compose affine maps: a -> eps x^s a(x^t); iterate on identity-ish: represent as (eps,s,t) composition
    cur=(eps,s,t); k=1
    while not (cur[0]==1 and cur[1]%N==0 and cur[2]%N==1):
        e1,s1,t1=cur; # compose cur after (eps,s,t): a -> eps x^s ( cur(a) )(x^t)... just iterate numerically
        # (T a)(x)=eps x^s a(x^t); T^2 a = eps x^s (eps (x^t)^s a((x^t)^t)) = x^{s+ts} a(x^{t^2})
        cur=(e1*eps, (s1 + t1*s)%N if False else (s + t*s1)%N, (t*t1)%N); k+=1
        if k>10000: return None
    return k
if __name__=="__main__":
    for v in (3,5,7,9,15,21,35):
        N=8*v; res={}
        for t in range(1,N,2):
            if gcd(t,N)!=1: continue
            for s in range(N):
                for eps in (1,-1):
                    c=count_inv(v,[(eps,s,t)])
                    if c>0 and (t%N!=1 or s!=0 or eps!=1):
                        key=(t%8,t%v if v>1 else 0)
                        res.setdefault(key,[]).append((eps,s,t,c))
        print('v=%d: elements (eps,s,t) != id admitting an invariant skew sequence, grouped by (t mod 8, t mod v):'%v)
        for key in sorted(res):
            lst=res[key]
            hh=key[1]
            o=pow(hh,1,v)
            # describe h order in Z_v^*
            ord_h=1; x=hh%v
            while x!=1%v and ord_h<1000: x=(x*hh)%v; ord_h+=1
            print('   t mod 8 = %d, h = t mod v = %s (ord %s, h==-1: %s): %d elements, e.g. %s'%(key[0],hh,ord_h,hh==v-1,len(lst),lst[:2]))

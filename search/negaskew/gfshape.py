"""Multiplier shapes over the Galois field GF(p^k) as the additive group G = F_q (q = p^k), the quaternary skew
Goethals-Seidel problem of README.md, H = a subgroup of F_q* of odd order acting by multiplication.  Writes the shape
file of (G, H) (the same format as qdrv.shape) and runs the exhaustive join `qskew2` on it (mode 1, first solution).

usage: python3 gfshape.py                       the six shapes of NOTES.md section 88 (361, 841, 1369, 1849, 2197 twice)
       python3 gfshape.py p k |H| [p k |H| ...]  other shapes, e.g.  python3 gfshape.py 19 2 45   (q = 361, 9 orbits)
Needs ./qskew2 (see README.md) and sympy; shape files go to shapes/gf_<p>_<k>_<|H|>.txt.  Expected: COUNT 0 on every shape."""
import numpy as np, itertools, sys, subprocess, os
def irreducible(p,k):
    # find monic irreducible poly of degree k over F_p by brute force (check no roots / factors via order test of x)
    from sympy import Poly, symbols, GF
    x=symbols('x')
    for coeffs in itertools.product(range(p),repeat=k):
        P=Poly([1]+list(coeffs),x,modulus=p)
        if P.is_irreducible: return [int(c) for c in reversed(coeffs)]  # low-degree-first for coefficients of x^0..x^{k-1}
def gf(p,k):
    mod=irreducible(p,k)  # x^k = -sum mod[i] x^i
    q=p**k
    def toidx(t):
        r=0
        for d in reversed(t): r=r*p+d
        return r
    def totup(i):
        t=[]
        for _ in range(k): t.append(i%p); i//=p
        return t
    def mul(a,b):
        A=totup(a);B=totup(b); prod=[0]*(2*k-1)
        for i in range(k):
            for j in range(k): prod[i+j]=(prod[i+j]+A[i]*B[j])%p
        for d in range(2*k-2,k-1,-1):
            c=prod[d]
            if c:
                prod[d]=0
                for i in range(k): prod[d-k+i]=(prod[d-k+i]-c*mod[i])%p
        return toidx(prod[:k])
    def add(a,b): return toidx([(x+y)%p for x,y in zip(totup(a),totup(b))])
    def neg(a): return toidx([(-x)%p for x in totup(a)])
    # primitive element
    for g in range(2,q):
        x=1;o=0
        while True:
            x=mul(x,g);o+=1
            if x==1: break
        if o==q-1: break
    return q,mul,add,neg,g
def gfshape(p,k,hord,fname):
    q,mul,add,neg,g=gf(p,k)
    e=(q-1)//hord
    h=1
    for _ in range(e): h=mul(h,g)
    # H = <h> of order hord ; orbits of F_q under multiplication by H
    Hs=[1]; x=h
    while x!=1: Hs.append(x); x=mul(x,h)
    assert len(Hs)==hord and (neg(1) not in Hs)
    orb=[-1]*q; reps=[]
    for a in range(q):
        if orb[a]<0:
            reps.append(a)
            for hh in Hs: orb[mul(a,hh)]=len(reps)-1
    n=len(reps); sz=[orb.count(i) for i in range(n)]; ng=[orb[neg(r)] for r in reps]
    N=np.zeros((n,n,n),dtype=np.int64)
    for c in range(1,n):
        s=reps[c]
        for a in range(q): N[c,orb[a],orb[add(a,s)]]+=1
    with open(fname,'w') as f:
        f.write(f"{q} {n}\n"+" ".join(map(str,sz))+"\n"+" ".join(map(str,ng))+"\n")
        for c in range(1,n):
            for a in range(n): f.write(" ".join(map(str,N[c,a]))+"\n")
    return q,n
if __name__=="__main__":
    args=[int(a) for a in sys.argv[1:]]
    if args: assert len(args)%3==0, "usage: gfshape.py [p k |H| ...]"
    cases=[tuple(args[i:i+3]) for i in range(0,len(args),3)] if args else [(37,2,171),(43,2,231),(13,3,549),(13,3,183),(19,2,45),(29,2,105)]
    os.makedirs('shapes',exist_ok=True)
    for (p,k,h) in cases:
        fn=f"shapes/gf_{p}_{k}_{h}.txt"; q,n=gfshape(p,k,h,fn)
        for mode in (1,):
            r=subprocess.run(['./qskew2',fn,str(mode),'1','1'],capture_output=True,text=True,timeout=int(os.environ.get('TMO','600')))
            print(q,'H',h,'n',n,'mode',mode,[l for l in r.stdout.split('\n') if l.startswith(('SOL','COUNT','INCOMPLETE'))],flush=True)

import numpy as np, subprocess, sys, json, itertools, os, math
from math import gcd
def subgroup(v,gens):
    H={1}
    frontier=[1]
    while frontier:
        x=frontier.pop()
        for g in gens:
            y=(x*g)%v
            if y not in H: H.add(y); frontier.append(y)
    return sorted(H)
def orbits(v,H):
    orb=[-1]*v; reps=[]
    # order orbits: orbit of 0 first, then by smallest element
    for g in range(v):
        if orb[g]<0:
            idx=len(reps); reps.append(g)
            for h in H: orb[(g*h)%v]=idx
    return orb,reps
def shape(v,gens,fname=None):
    H=subgroup(v,gens)
    assert (v-1) not in H or v==1, "-1 in H"
    orb,reps=orbits(v,H); n=len(reps)
    sz=[orb.count(i) for i in range(n)]
    neg=[orb[(-reps[i])%v] for i in range(n)]
    # N tensors
    N=np.zeros((n,n,n),dtype=np.int64)
    for c in range(1,n):
        s=reps[c]
        for g in range(v):
            N[c,orb[g],orb[(g+s)%v]]+=1
    d=dict(v=v,n=n,H=H,orb=orb,reps=reps,sz=sz,neg=neg,N=N)
    if fname:
        with open(fname,'w') as f:
            f.write(f"{v} {n}\n"+" ".join(map(str,sz))+"\n"+" ".join(map(str,neg))+"\n")
            for c in range(1,n):
                for a in range(n): f.write(" ".join(map(str,N[c,a]))+"\n")
    return d
def run(v,gens,mode,maxprint=3,firstonly=0,tag='t',timeout=None):
    fn=f"shapes/{tag}_{v}.txt"; os.makedirs('shapes',exist_ok=True)
    d=shape(v,gens,fn)
    cmd=['./qskew',fn,str(mode),str(maxprint),str(firstonly)]
    r=subprocess.run(cmd,capture_output=True,text=True,timeout=timeout)
    return d,r.stdout,r.stderr
def expand(d,e):
    return [int(e[d['orb'][g]]) for g in range(d['v'])]
# ---- real 8v matrix from quaternary exponents (G level, Z_v) ----
def real_matrix(v,E):  # E: list of 4 exponent lists (length v) for A,B,C,D (units i^e)
    ab=[(1,1),(-1,1),(-1,-1),(1,-1)]  # (1+i)*i^e = a+ib
    def L(e):
        M=np.zeros((2*v,2*v),dtype=np.float32)
        for g in range(v):
            for h in range(v):
                a,b=ab[e[(h-g)%v]]
                M[2*g:2*g+2,2*h:2*h+2]=[[a,-b],[b,a]]
        return M
    Ls=[L(e) for e in E]
    R=np.zeros((2*v,2*v),dtype=np.float32)
    for g in range(v):
        R[2*g,2*((-g)%v)]=1; R[2*g+1,2*((-g)%v)+1]=-1
    A,B,C,D=Ls
    T=lambda X:X.T
    rows=[[A,B@R,C@R,D@R],[-B@R,A,T(D)@R,-T(C)@R],[-C@R,-T(D)@R,A,T(B)@R],[-D@R,T(C)@R,-T(B)@R,A]]
    H=np.block(rows)
    return H
def check_real(H):
    n=H.shape[0]
    ok_pm=bool(np.all(np.abs(H)==1))
    P=H@H.T
    return ok_pm, bool(np.array_equal(P,n*np.eye(n,dtype=np.float32))), bool(np.array_equal(H+H.T,2*np.eye(n,dtype=np.float32)))
if __name__=="__main__":
    v=int(sys.argv[1]); mode=int(sys.argv[2]); gens=json.loads(sys.argv[3]) if len(sys.argv)>3 else []
    d,out,err=run(v,gens,mode,maxprint=3)
    print(err.strip()); print(out.strip())
    for line in out.split('\n'):
        if line.startswith('SOL'):
            parts=line.split()[1:]
            E=[expand(d,[int(c) for c in p]) for p in parts]
            if mode==0: E=[[2*(x//2) if False else x for x in e] for e in E]
            H=real_matrix(v,E) if mode==1 else None
            if H is not None: print('real order',H.shape[0],check_real(H))

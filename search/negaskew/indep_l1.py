"""Independent level-1 join of the quaternary Paley plug-in (written from the Fourier derivation, sharing no code with qp.c): the check of the engine's
level 1 (NOTES.md section 89).  G = Z_m x F_q (q = 3 mod 4 prime), blocks constant on (K-orbit of H_K = <hK>, {0 | squares | non-squares}) with unit values;
level 1:  w(k) = x0(k) + (q-1)(x1(k)+x2(k))/2,  sum_X PAF_K(w_X) = 4 m q delta_K.  Blocks B, C, D: x0(0) = 1 (unit normalisation); block A twisted skew
u(-g) = i conj u(g) (g != 0), u(0) in {1, -i}.  The only necessary filter is the weight W_X = sum_k c_X(k) <= 4m (|w|^2 = 1 + (q-1) c); --psd adds the
single-block bound |what|^2 <= 4 m q at every character.  Lists, hash join on the exact PAF vectors, exact verification of every match.
usage: python3 indep_l1.py m q "[H_K generators]" [--psd]       numpy only; writes runs/indep_l1_<m>_<q>_<|H_K|>.json (the instances)
Example: python3 indep_l1.py 13 23 "[3]" --psd  ->  LEVEL1_INSTANCES 36000 (the 36,000 level-1 instances of qp at v = 299, 5 s);  19 31 "[7]" --psd gives the 261,984 of v = 589 (13 minutes)."""
import sys, itertools, json, time
import numpy as np
np.seterr(over='ignore')     # the hash sums of uint64 wrap around on purpose
m=int(sys.argv[1]); q=int(sys.argv[2]); hK=json.loads(sys.argv[3]); usepsd=('--psd' in sys.argv)
U=[1,1j,-1,-1j]
# K orbits
H={1}; fr=[1]
while fr:
    x=fr.pop()
    for g in hK:
        y=(x*g)%m
        if y not in H: H.add(y); fr.append(y)
assert (m-1) not in H
orb=[-1]*m; reps=[]
for k in range(m):
    if orb[k]<0:
        reps.append(k)
        for h in H: orb[(k*h)%m]=len(reps)-1
nu=len(reps); ell=[orb.count(o) for o in range(nu)]; neg=[orb[(-reps[o])%m] for o in range(nu)]
print('m',m,'q',q,'|H_K|',len(H),'nu',nu,'ell',ell,'neg',neg,flush=True)
# symbols: distinct (x0, esum=x1+x2): store x0 (unit exp), and E2 = x1+x2 as complex int, 2w = 2*x0+(q-1)*E2
syms={}
for a in range(4):
    for b in range(4):
        for c in range(4):
            E2=U[b]+U[c]
            syms[(a,complex(E2))]=1
syms=sorted(syms.keys(),key=lambda t:(t[0],t[1].real,t[1].imag))
print('number of (x0,x1+x2) classes',len(syms))
def w2(s): a,E2=s; return 2*U[a]+(q-1)*E2
def cweight(s):
    W=w2(s); n2=(W.real**2+W.imag**2)/4.0   # |w|^2
    c=(n2-1)/(q-1); assert abs(c-round(c))<1e-9 and round(c)>=0,(s,c); return int(round(c))
S=[(w2(s),cweight(s),s) for s in syms]
cap=4*m
# per-orbit options
def gen_free():
    out=[]
    def rec(o,Wt,cur):
        if o==nu: out.append(list(cur)); return
        for (w,c,s) in S:
            if o==0 and s[0]!=0: continue
            W2=Wt+ell[o]*c
            if W2>cap: continue
            cur.append(w); rec(o+1,W2,cur); cur.pop()
    rec(0,0,[])
    return out
def gen_skew():
    out=[]
    # origin symbols: x0 in {1,-i} (exp 0 or 3), x1+x2 = +-(1+i)  (x2 = i conj x1)
    org=[s for s in S if s[2][0] in (0,3) and s[2][1] in (1+1j,-1-1j)]
    # check origin legality directly: exists (x1,x2) units with x2 = i*conj(x1) and x1+x2=E2
    def legal(E2):
        return any(abs(U[c]-1j*np.conj(U[b]))<1e-9 and abs(U[b]+U[c]-E2)<1e-9 for b in range(4) for c in range(4))
    org=[s for s in org if legal(s[2][1])]
    # but also all origin symbols satisfying skew (not only those two E2) -- enumerate by skew relation to be safe
    org=[]
    for a in (0,3):
        for b in range(4):
            c=[t for t in range(4) if abs(U[t]-1j*np.conj(U[b]))<1e-9][0]
            E2=U[b]+U[c]
            for s in S:
                if s[2][0]==a and abs(s[2][1]-E2)<1e-9: org.append(s)
    org=list({id(s):s for s in org}.values())
    pairs=[o for o in range(1,nu) if o<neg[o]]
    def rec(i,Wt,cur,base):
        if i==len(pairs):
            out.append(base); 
            return
    res=[]
    def rec2(i,Wt,asg,w0):
        if i==len(pairs):
            res.append((w0,dict(asg))); return
        o=pairs[i]
        for (w,c,s) in S:
            W2=Wt+2*ell[o]*c
            if W2>cap: continue
            asg[o]=w; rec2(i+1,W2,asg,w0); 
        asg.pop(o,None) if False else None
    for (w0,c0,s0) in org:
        if c0>cap: continue
        rec2(0,c0,{},w0)
    # build full 2w vectors per orbit
    full=[]
    for (w0,asg) in res:
        v=[0]*nu; v[0]=w0
        for o,w in asg.items():
            v[o]=w; v[neg[o]]=1j*np.conj(w)
        full.append(v)
    return full
LA=gen_skew(); LF=gen_free()
print('|A|=',len(LA),'|L|=',len(LF),flush=True)
# PAF vectors (of 2w) over the nonzero lag class reps and lag 0
lagk=[reps[o] for o in range(1,nu)]
kidx=np.array(orb)
def paf_vec(v):
    w=np.array(v,dtype=complex)[kidx]    # length m, 2w
    out=[np.sum(np.abs(w)**2)]
    for s in lagk:
        out.append(np.sum(w*np.conj(w[(np.arange(m)+s)%m])))
    return out
def vecs(L):
    X=np.zeros((len(L),1+2*(nu-1)),dtype=np.int64)
    for i,v in enumerate(L):
        p=paf_vec(v)
        X[i,0]=int(round(p[0].real))
        for j in range(nu-1):
            assert abs(p[1+j].real-round(p[1+j].real))<1e-6 and abs(p[1+j].imag-round(p[1+j].imag))<1e-6
            X[i,1+2*j]=int(round(p[1+j].real)); X[i,2+2*j]=int(round(p[1+j].imag))
    return X
if usepsd:
    lim=4*m*q
    def psdok(v):
        w=np.array(v,dtype=complex)[kidx]/2
        fw=np.fft.fft(w); return np.all(np.abs(fw)**2<=lim+1e-6)
    LA=[v for v in LA if psdok(v)]; LF=[v for v in LF if psdok(v)]
    print('after psd: |A|=',len(LA),'|L|=',len(LF))
XA=vecs(LA); XF=vecs(LF)
T=np.zeros(1+2*(nu-1),dtype=np.int64); T[0]=16*m*q   # sum of PAF(2w)(0)= 4*4mq
rng=np.random.default_rng(12345)
coef=rng.integers(1,2**62,size=XA.shape[1],dtype=np.int64).astype(np.uint64)|np.uint64(1)
def hsh(X): return (X.astype(np.uint64)*coef[None,:]).sum(axis=1,dtype=np.uint64)
hA=hsh(XA); hF=hsh(XF); hT=hsh(T[None,:])[0]
t0=time.time()
# left = A + B(F), right = C + D
left=(hA[:,None]+hF[None,:]).ravel()
order=np.argsort(left,kind='stable'); ls=left[order]
cnt=0; matches=[]
nF=len(LF)
for c in range(nF):
    need=(hT-hF[c]-hF).astype(np.uint64)         # needed left key for each d
    pos=np.searchsorted(ls,need,side='left'); pos2=np.searchsorted(ls,need,side='right')
    mk=pos2>pos
    for d in np.nonzero(mk)[0]:
        for p in range(pos[d],pos2[d]):
            idx=order[p]; a=idx//nF; b=idx%nF
            matches.append((a,b,c,int(d)))
print('hash matches',len(matches),'time',round(time.time()-t0,1),flush=True)
# exact verification of every match
good=0
for (a,b,c,d) in matches:
    tot=XA[a]+XF[b]+XF[c]+XF[d]
    if np.array_equal(tot,T): good+=1
print('LEVEL1_INSTANCES (ordered, B,C,D normalised x0(0)=1):',good,' (hash false positives:',len(matches)-good,')')
import os
os.makedirs('runs',exist_ok=True)
json.dump([[int(x) for x in t] for t in matches if np.array_equal(XA[t[0]]+XF[t[1]]+XF[t[2]]+XF[t[3]],T)],open(f'runs/indep_l1_{m}_{q}_{len(H)}.json','w'))

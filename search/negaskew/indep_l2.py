"""Independent level 1 + level 2 (the full completion count) of the quaternary Paley plug-in, from the Fourier derivation (sharing no code with qp.c): the check
of the engine's both levels (NOTES.md section 89).  Level 2 (the characters a != 0 of F_q) with s = x0 - e, t = (x1-x2)/2, e = (x1+x2)/2:
  sum_X [ PAF_K(s_X) + q PAF_K(t_X) ] = 4 m q delta,   sum_X [ C_{t,s} - C_{s,t} ] = 0     (K-correlations, exact Gaussian integers).
usage: python3 indep_l2.py m q "[H_K generators]" [--sample N seed] [--out sols.txt]      numpy only; --sample: level 2 on N random level-1 instances
Prints the level-1 instances and 'LEVEL2: instances processed .. completing instances .. solutions ..'; --out writes the solutions as SOL lines in the format of qp
(four blocks x0|x1|x2).  Example: python3 indep_l2.py 13 23 "[3]" --out runs/p299_indep.txt  ->  36,000 level-1 instances, 4,800 completing, 13,056 solutions (6 minutes),
the same set of solutions as  ./qp shapes/p299.txt l2c=1 maxprint=20000;  python3 indep_l2.py 13 23 "[3]" --sample 200 1  is the quick version."""
import sys, itertools, json, time, random
import numpy as np
np.seterr(over='ignore')     # the hash sums of uint64 wrap around on purpose
m=int(sys.argv[1]); q=int(sys.argv[2]); hK=json.loads(sys.argv[3])
sample=None; seed=1; outf=None
if '--sample' in sys.argv: i=sys.argv.index('--sample'); sample=int(sys.argv[i+1]); seed=int(sys.argv[i+2])
if '--out' in sys.argv: outf=sys.argv[sys.argv.index('--out')+1]
U=[1,1j,-1,-1j]
H={1}; fr=[1]
while fr:
    x=fr.pop()
    for g in hK:
        y=(x*g)%m
        if y not in H: H.add(y); fr.append(y)
orb=[-1]*m; reps=[]
for k in range(m):
    if orb[k]<0:
        reps.append(k)
        for h in H: orb[(k*h)%m]=len(reps)-1
nu=len(reps); ell=[orb.count(o) for o in range(nu)]; neg=[orb[(-reps[o])%m] for o in range(nu)]
kidx=np.array(orb)
# symbols (x0 exponent a, E2 = x1+x2)
SY=sorted({(a,complex(U[b]+U[c])) for a in range(4) for b in range(4) for c in range(4)},key=lambda t:(t[0],t[1].real,t[1].imag))
def w2(s): a,E2=s; return 2*U[a]+(q-1)*E2
def cw(s):
    W=w2(s); c=(abs(W)**2/4-1)/(q-1); assert abs(c-round(c))<1e-9 and round(c)>=0; return int(round(c))
SW=[(w2(s),cw(s)) for s in SY]
cap=4*m
def conjE(E2): return 1j*np.conj(E2)
def partner_sym(s):  # (x0,E2) -> (i conj x0, i conj E2)
    a,E2=s; a2=[t for t in range(4) if abs(U[t]-1j*np.conj(U[a]))<1e-9][0]
    E2p=1j*np.conj(E2); return SY.index((a2,complex(round(E2p.real),round(E2p.imag))))
def gen_free():
    out=[]
    def rec(o,Wt,cur):
        if o==nu: out.append(tuple(cur)); return
        for si,(s) in enumerate(SY):
            if o==0 and s[0]!=0: continue
            W2=Wt+ell[o]*SW[si][1]
            if W2>cap: continue
            cur.append(si); rec(o+1,W2,cur); cur.pop()
    rec(0,0,[]); return out
def gen_skew():
    org=[]
    for a in (0,3):
        for b in range(4):
            c=[t for t in range(4) if abs(U[t]-1j*np.conj(U[b])*1j)<1e-9]
            # x2 = i conj(x1)
            c=[t for t in range(4) if abs(U[t]-1j*np.conj(U[b]))<1e-9][0]
            E2=U[b]+U[c]
            si=SY.index((a,complex(round(E2.real),round(E2.imag)))); 
            if si not in org: org.append(si)
    pairs=[o for o in range(1,nu) if o<neg[o]]
    res=[]
    def rec(i,Wt,asg,s0):
        if i==len(pairs): res.append((s0,dict(asg))); return
        o=pairs[i]
        for si in range(len(SY)):
            W2=Wt+2*ell[o]*SW[si][1]
            if W2>cap: continue
            asg[o]=si; rec(i+1,W2,asg,s0)
        asg.pop(o,None)
    for s0 in org:
        if SW[s0][1]>cap: continue
        rec(0,SW[s0][1],{},s0)
    full=[]
    for (s0,asg) in res:
        v=[0]*nu; v[0]=s0
        for o,si in asg.items(): v[o]=si; v[neg[o]]=partner_sym(SY[si])
        full.append(tuple(v))
    return full
LA=gen_skew(); LF=gen_free()
lagk=[reps[o] for o in range(1,nu)]
def paf_vec(sy):
    w=np.array([w2(SY[si]) for si in sy],dtype=complex)[kidx]
    out=[np.sum(np.abs(w)**2)]
    for s in lagk: out.append(np.sum(w*np.conj(w[(np.arange(m)+s)%m])))
    return out
def vecs(L):
    X=np.zeros((len(L),1+2*(nu-1)),dtype=np.int64)
    for i,sy in enumerate(L):
        p=paf_vec(sy); X[i,0]=int(round(p[0].real))
        for j in range(nu-1): X[i,1+2*j]=int(round(p[1+j].real)); X[i,2+2*j]=int(round(p[1+j].imag))
    return X
lim=4*m*q
def psdok(sy):
    w=np.array([w2(SY[si]) for si in sy],dtype=complex)[kidx]/2
    return bool(np.all(np.abs(np.fft.fft(w))**2<=lim+1e-6))
LA=[x for x in LA if psdok(x)]; LF=[x for x in LF if psdok(x)]
print('|A|',len(LA),'|L|',len(LF),flush=True)
XA=vecs(LA); XF=vecs(LF)
T=np.zeros(XA.shape[1],dtype=np.int64); T[0]=16*m*q
rng=np.random.default_rng(12345); coef=rng.integers(1,2**62,size=XA.shape[1],dtype=np.int64).astype(np.uint64)|np.uint64(1)
def hsh(X): return (X.astype(np.uint64)*coef[None,:]).sum(axis=1,dtype=np.uint64)
hA=hsh(XA); hF=hsh(XF); hT=hsh(T[None,:])[0]
left=(hA[:,None]+hF[None,:]).ravel(); order=np.argsort(left,kind='stable'); ls=left[order]; nF=len(LF)
inst=[]
for c in range(nF):
    need=(hT-hF[c]-hF).astype(np.uint64)
    p1=np.searchsorted(ls,need,'left'); p2=np.searchsorted(ls,need,'right')
    for d in np.nonzero(p2>p1)[0]:
        for p in range(p1[d],p2[d]):
            idx=order[p]; a=idx//nF; b=idx%nF
            if np.array_equal(XA[a]+XF[b]+XF[c]+XF[d],T): inst.append((a,b,c,int(d)))
print('level-1 instances',len(inst),flush=True)
# ---------------- level 2
def unit_idx(z): return [t for t in range(4) if abs(U[t]-z)<1e-9][0]
def options(E2,skew_origin=False):
    out=[]
    for b in range(4):
        for c in range(4):
            if abs(U[b]+U[c]-E2)<1e-9:
                if skew_origin and abs(U[c]-1j*np.conj(U[b]))>1e-9: continue
                out.append((b,c))
    return out
def fibre(sy,kind):
    """all (x1,x2) exponent assignments for the orbits that carry free choices; returns arrays s2 (m), list of t2 (n_el, m) and the digit table"""
    sval=[2*U[SY[si][0]]-SY[si][1] for si in sy]            # s2 per orbit
    free=[o for o in range(nu) if (kind==1 or o==0 or o<neg[o])]
    opts=[options(SY[sy[o]][1],skew_origin=(kind==0 and o==0)) for o in free]
    combos=list(itertools.product(*opts))
    n_el=len(combos)
    X1=np.zeros((n_el,nu),dtype=np.int64); X2=np.zeros((n_el,nu),dtype=np.int64)
    for j,o in enumerate(free):
        X1[:,o]=[cmb[j][0] for cmb in combos]; X2[:,o]=[cmb[j][1] for cmb in combos]
    if kind==0:
        for o in range(1,nu):
            if o<neg[o]:
                p=neg[o]    # partner: x1' = i conj x2, x2' = i conj x1  (exponent of i conj(u^e) = 1 - e)
                X1[:,p]=(1-X2[:,o])%4; X2[:,p]=(1-X1[:,o])%4
    Ur=np.array(U)
    t2=(Ur[X1]-Ur[X2])           # (n_el, nu) = x1 - x2  (=2t)
    return np.array(sval,dtype=complex),t2,X1,X2
lagsK=[0]+lagk
def block_vec(sv,t2):
    s=sv[kidx]; t=t2[:,kidx]                  # s: (m), t: (n,m)
    n_el=t.shape[0]; out=np.zeros((n_el,1+8+1+8),dtype=np.int64)
    for j,sh in enumerate(lagsK):
        ix=(np.arange(m)+sh)%m
        ss=np.sum(s*np.conj(s[ix])); tt=np.sum(t*np.conj(t[:,ix]),axis=1)
        d1=ss+q*tt                                   # PAF(s2)+q PAF(t2)  (= 4 x true)
        cts=np.sum(t*np.conj(s[ix]),axis=1); cst=np.sum(s*np.conj(t[:,ix]),axis=1)
        d2=cts-cst
        if j==0: out[:,0]=np.rint(d1.real).astype(np.int64); out[:,9]=np.rint(d2.imag).astype(np.int64)
        else:
            out[:,1+2*(j-1)]=np.rint(d1.real); out[:,2+2*(j-1)]=np.rint(d1.imag)
            out[:,10+2*(j-1)]=np.rint(d2.real); out[:,11+2*(j-1)]=np.rint(d2.imag)
    return out
coef2=rng.integers(1,2**62,size=18,dtype=np.int64).astype(np.uint64)|np.uint64(1)
def hsh2(X): return (X.astype(np.uint64)*coef2[None,:18]).sum(axis=1,dtype=np.uint64)
T2=np.zeros(18,dtype=np.int64); T2[0]=16*m*q
hT2=hsh2(T2[None,:])[0]
def strs(sy,X1,X2,i):
    s0=''.join(str(SY[si][0]) for si in sy); s1=''.join(str(int(X1[i,o])) for o in range(nu)); s2=''.join(str(int(X2[i,o])) for o in range(nu))
    return f"{s0}|{s1}|{s2}"
sel=list(range(len(inst)))
if sample: random.Random(seed).shuffle(sel); sel=sorted(sel[:sample])
sols=[]; nhit=0; t0=time.time(); cache={}
def get_block(kind,idx):
    key=(kind,idx)
    if key in cache: return cache[key]
    sy=(LA if kind==0 else LF)[idx]
    sv,t2,X1,X2=fibre(sy,kind); V=block_vec(sv,t2)
    r=(sy,X1,X2,V,hsh2(V)); 
    if len(cache)>4000: cache.clear()
    cache[key]=r; return r
for n_i,ii in enumerate(sel):
    a,b,c,d=inst[ii]
    BA=get_block(0,a); BB=get_block(1,b); BC=get_block(1,c); BD=get_block(1,d)
    L=(BA[4][:,None]+BB[4][None,:]).ravel(); R=(BC[4][:,None]+BD[4][None,:]).ravel()
    # need L + R = hT2
    need=(hT2-R).astype(np.uint64)
    Ls=np.sort(L)
    pos=np.searchsorted(Ls,need); pos[pos>=len(Ls)]=len(Ls)-1
    cand=np.nonzero(Ls[pos]==need)[0]
    found=0
    if len(cand):
        # exact verification of every candidate
        Lidx={}
        for i,hv in enumerate(L): Lidx.setdefault(int(hv),[]).append(i)
        for r in cand:
            for l in Lidx[int(need[r])]:
                i1,i2=divmod(l,len(BB[4])); i3,i4=divmod(int(r),len(BD[4]))
                tot=BA[3][i1]+BB[3][i2]+BC[3][i3]+BD[3][i4]
                if np.array_equal(tot,T2):
                    found+=1; sols.append(' '.join(strs(B[0],B[1],B[2],i) for B,i in ((BA,i1),(BB,i2),(BC,i3),(BD,i4))))
    nhit+= (found>0)
    if n_i%500==499: print('instances done',n_i+1,'of',len(sel),'solutions so far',len(sols),'time',round(time.time()-t0),flush=True)
print('LEVEL2: instances processed',len(sel),'completing instances',nhit,'solutions',len(sols),'time',round(time.time()-t0,1))
if outf: open(outf,'w').write('\n'.join('SOL '+x for x in sols)+'\n')

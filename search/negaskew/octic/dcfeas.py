"""DC (row-sum) feasibility of octic shapes: exists R_A (skew type), R_B, R_C, R_D with sum_X |R_X|^2 = 16 v in Z[sqrt2].
|R|^2 = P + sqrt2 Q, P = sum r_i^2, Q = r0r1+r1r2+r2r3-r0r3.  R_B,R_C,R_D: coordinates independent subset sums of orbit sizes; R_A = (1,r1,r2,r1).
A block constant on the H-orbits (|H| = h acting freely off 0) has R = U(0) + h S with S in Z[w] of even coordinates, so sum_X |R_X|^2 = 16 v needs
(2h - 1)^2 + 15 <= 16 v: h <= 2 sqrt(v) + 1/2 (NOTES.md section 89).
usage: python3 targets.py && python3 shapes_oct.py NMAX && python3 dcfeas.py NMAX      NMAX = largest number of orbits (65 for the full table)
 reads runs/oct_shapes_n<NMAX>.json (shapes_oct.py), writes runs/dcfeas_n<NMAX>.json and prints one line per shape: dc_infeasible,
 dc_feasible, dead_selfneg (an orbit with -O = O cannot carry a skew block) or dense_not_evaluated.  Every shape with at most 14 orbits of
 the 61 targets is infeasible; the feasible ones (n = 17 at 337, 21 at 101 and 821, 23 at 947, ...) need lists of at least 2^57 blocks.
module: orbit_data(v, gens), subset_sums, feas(v, sz, neg)"""
import numpy as np, json, sys, itertools, os
from math import gcd, isqrt
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from shapes_oct import gen_sub
def orbit_data(v,gens):
    H=sorted(gen_sub(v,gens)); seen=[-1]*v; reps=[]
    for g in range(v):
        if seen[g]<0:
            i=len(reps); reps.append(g)
            for h in H: seen[g*h%v]=i
    n=len(reps); sz=[seen.count(i) for i in range(n)]; neg=[seen[(-reps[i])%v] for i in range(n)]
    return H,seen,reps,sz,neg
def subset_sums(sizes,start=(1,-1),M=None):
    M=M or sum(sizes)+5
    cur=np.zeros(2*M+1,bool)
    for s in start: cur[M+s]=True
    for s in sizes:
        new=np.zeros_like(cur)
        new[s:]|=cur[:-s] if s>0 else cur
        new[:-s]|=cur[s:] if s>0 else cur
        cur=new
    return np.nonzero(cur)[0]-M
def feas(v,sz,neg,maxSL=4000,verbose=False):
    n=len(sz); T=16*v; B=isqrt(T)+1
    if any(neg[a]==a for a in range(1,n)): return 'dead_selfneg',None
    Z=subset_sums(sz[1:],start=(1,-1),M=v+5)       # orbit 0 contributes +-1
    Zc=Z[np.abs(Z)<=B]
    if len(Zc)==0: return 'dc_infeasible',{'SL':0}
    if len(Zc)**4>2.5e7: return 'dense_not_evaluated',{'Z':len(Zc)}
    r=np.array(np.meshgrid(Zc,Zc,Zc,Zc,indexing='ij')).reshape(4,-1).astype(np.int64)
    P=(r**2).sum(0); Q=r[0]*r[1]+r[1]*r[2]+r[2]*r[3]-r[0]*r[3]
    ok=(P<=T)&((P*P-2*Q*Q)>=0)&(P+np.sqrt(2)*Q>=-1e-9)&(P-np.sqrt(2)*Q>=-1e-9)
    P=P[ok];Q=Q[ok]; K=int(4*T+10)
    enc=np.unique(P*(2*K)+(Q+K))
    # A-type
    pairs=[a for a in range(1,n) if a<neg[a]]; ps=[sz[a] for a in pairs]
    Z1=subset_sums([2*s for s in ps]+[],start=(1,-1),M=v+5) if ps else np.array([1,-1])
    # r1 = c + sum sz*t, t in {-2,0,2}: sums of {-2s,0,2s}: use subset sums of {+-s}*... handle: t=c1+c3 in {-2,0,2}
    R1=set([1,-1])
    for s in ps: R1={x+t*s for x in R1 for t in (-2,0,2)}
    R2=set([1,-1])
    for s in ps: R2={x+t*s for x in R2 for t in (-2,2)}
    SA=set()
    for a in R1:
        for b in R2:
            Pa=1+2*a*a+b*b; Qa=2*a*b
            if Pa<=T and Pa+np.sqrt(2)*Qa>=-1e-9 and Pa-np.sqrt(2)*Qa>=-1e-9: SA.add((Pa,Qa))
    if not SA: return 'dc_infeasible',{'SA':0,'SL':len(enc)}
    if len(enc)>maxSL: return 'dense_not_evaluated',{'SL':len(enc),'SA':len(SA)}
    PL=enc//(2*K); QL=enc%(2*K)-K
    # SLL = SL+SL
    PP=(PL[:,None]+PL[None,:]).ravel(); QQ=(QL[:,None]+QL[None,:]).ravel()
    m=PP<=T; sll=np.unique(PP[m]*(2*K)+(QQ[m]+K))
    # note Q-sums can exceed K range; mask by |Q|<K
    for (Pa,Qa) in SA:
        Pr=T-Pa; Qr=-Qa
        # need b + (c+d) = (Pr,Qr): b in SL
        Pc=Pr-PL; Qc=Qr-QL
        m2=(Pc>=0)&(np.abs(Qc)<K)
        if m2.any():
            e=Pc[m2]*(2*K)+(Qc[m2]+K)
            if np.isin(e,sll).any(): return 'dc_feasible',{'SL':len(enc),'SA':len(SA)}
    return 'dc_infeasible',{'SL':len(enc),'SA':len(SA)}
if __name__=="__main__":
    nmax=int(sys.argv[1]); d=json.load(open(f'runs/oct_shapes_n{nmax}.json')); out=[]
    for v,(m,sh) in sorted(d.items(),key=lambda x:int(x[0])):
        v=int(v)
        for n,o,g in sh:
            H,seen,reps,sz,neg=orbit_data(v,g)
            st,info=feas(v,sz,neg)
            out.append(dict(v=v,m=m,n=n,order=o,gens=g,status=st,info=info))
            print(v,m,'n=%d h=%d'%(n,o),st,info,flush=True)
    json.dump(out,open(f'runs/dcfeas_n{nmax}.json','w'))

"""The border of the n = 955 core (k = 8, character set T = all 7 nontrivial characters of Q = Z2^3, pattern 2).

usage: python3 border8.py [--out DIR]

Border of complex-hadamard/TBORDERED.md for t = 7, q = 8: a 14 x 16 matrix Rhat and a 16 x 14 matrix Chat of quaternary entries whose
rows (columns) lie in V_T (each half is a zero-sum vector on Q = Z2^3: the Fourier transform vanishes off T) and are mutually
orthogonal, with Xb = -Rhat Mc* Chat / 16 quaternary or zero (a quaternary conference matrix of order 14). Mc is the 16 x 16 core on Q
built from the compressions of pattern 2, yA = yB = (0, 3), yC = yD = (-1, 1) (y(0), y(a) for a != 0), Mc = [[X, Y], [Y*, -X*]] with
X = alpha yA + conj(alpha) yB, Y = alpha yC + conj(alpha) yD, alpha = (1 + i)/2 (R = identity on Z2^3).
Ansatz: tau = the Singer cycle of GF(2^3) (companion matrix of x^3 + x + 1) acts on Q; the 7 rows of Rhat (and the 7 columns of Chat)
of one orbit are the tau-translates of one seed vector (the first 7 rows come from seed r1, the last 7 from r2; the same for columns).

Stages (run in this order, in memory; --out receives good_p1.npy, good_p2.npy, borders_p2.npy and border_k8_pattern2.json):
  enum     the 4900 zero-sum quaternary vectors of length 8 (exponents of i; as many 0s as 2s and as many 1s as 3s)
  sig      signature of a vector h = (<h, tau^j h>)_{j = 1,2,3}; seeds r = (h0 | h1) of 16 entries with h0[0] = 1 and
           signature(h1) = -signature(h0) (zero total autocorrelation for j = 1..6)
  good     pairs (r, c) of seeds with -r Mc* tau^k c / 16 zero or quaternary for all k = 0..6, for both patterns tried
           (pattern 1 = yA (0,1), yB (0,3), yC (1,1), yD (-1,3): no pair at all; pattern 2: 8820 pairs); writes good_p1.npy, good_p2.npy
  orth     orthogonality graphs on the row seeds and on the column seeds that occur in good pairs
  borders  all (r1, r2, c1, c2) with r1 perp r2 (all tau-translates), c1 perp c2, the zero pattern of Xb a permutation matrix:
           writes borders_p2.npy (28,812 of them)
  check    the first border: Rhat Rhat* = 16 I, Chat* Chat = 16 I, zero-sum halves, entries, Xb Xb* = 13 I
  save     border_k8_pattern2.json for the first border in the order of 'borders' (the file shipped in data/)
"""
import argparse, itertools, json, os, sys, time
from collections import Counter, defaultdict
import numpy as np
import tbio

ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("--out", default=None, help="output directory (default: current directory)")
args = ap.parse_args()
out = tbio.outdir(args.out)
T0 = time.time()
def say(*a):
    print("[%6.1fs]" % (time.time() - T0), *a); sys.stdout.flush()

# ---- Q = Z2^3 with the Singer cycle tau (companion of x^3 + x + 1); index a = 4*a0 + 2*a1 + a2
tau = np.array([[0,0,1],[1,0,1],[0,1,0]])
Qel = [np.array([(a>>2)&1,(a>>1)&1,a&1]) for a in range(8)]
def idx(v): return 4*int(v[0]%2)+2*int(v[1]%2)+int(v[2]%2)
pt = [idx(tau@Qel[a]) for a in range(8)]          # a -> tau a
perm = [list(range(8))]
for j in range(1,7):
    perm.append([pt[x] for x in perm[-1]])
perm16 = [perm[j]+[8+x for x in perm[j]] for j in range(7)]
neg = lambda s: tuple((-a,-b) for a,b in s)

# ---- stage enum: zero-sum quaternary 8-vectors
halves=[]
for e in itertools.product(range(4),repeat=8):
    c=[0]*4
    for x in e: c[x]+=1
    if c[0]==c[2] and c[1]==c[3]:
        halves.append(e)
H = np.array(halves,dtype=np.int8)      # 4900 x 8 exponents
Z = (1j)**H.astype(int)                  # complex
say("zero-sum quaternary vectors of length 8:", len(H))
assert len(H) == 4900

# ---- stage sig: signatures and seeds
sig = []
for h in range(len(H)):
    s=[]
    for j in (1,2,3):
        v = np.sum(Z[h]*np.conj(Z[h][perm[j]]))
        s.append((int(round(v.real)),int(round(v.imag))))
    sig.append(tuple(s))
groups=defaultdict(list)
for h,s in enumerate(sig): groups[s].append(h)
pairs=0
for s,l in groups.items():
    if neg(s) in groups: pairs+=len(l)*len(groups[neg(s)])
say("signatures",len(groups),"half-pairs (h0,h1) with zero total autocorrelation:",pairs)
seeds=[]
for s,l in groups.items():
    if neg(s) in groups:
        for h0 in l:
            if H[h0][0]!=0: continue
            for h1 in groups[neg(s)]:
                seeds.append((h0,h1))
seeds=np.array(seeds)
S = np.concatenate([Z[seeds[:,0]], Z[seeds[:,1]]],axis=1)   # n x 16 complex
say("seeds",S.shape)

# ---- the core on Q for the two compression patterns
al = (1+1j)/2
def devel(y):  # Q-developed matrix X[a,b] = y(a xor b); y given as (y0,u): value y0 at 0 else u
    M = np.zeros((8,8))
    for a in range(8):
        for b in range(8):
            M[a,b] = y[0] if a==b else y[1]
    return M
def Mc_from(y):
    yA,yB,yC,yD = [devel(v) for v in y]
    X = al*yA + np.conj(al)*yB
    Y = al*yC + np.conj(al)*yD      # R = identity on Z2^3
    return np.block([[X,Y],[Y.conj().T, -X.conj().T]])
patterns = {"p1": [(0,1),(0,3),(1,1),(-1,3)], "p2": [(0,3),(0,3),(-1,1),(-1,1)]}
Mc = Mc_from(patterns["p2"])
def quatok(v):
    """|v| = 0, or v in {1, i, -1, -i} (up to rounding)"""
    return (np.abs(v)<1e-6) | ((np.abs(np.abs(v)-1)<1e-6) & ((np.abs(v.real)<1e-6)|(np.abs(v.imag)<1e-6)))

# ---- stage good
good = {}
for name,y in patterns.items():
    Mcn = Mc_from(y)
    G = Mcn.conj().T@Mcn
    P_T = np.eye(16) - np.kron(np.eye(2), np.ones((8,8))/8)
    say(name, "Mc*Mc P_T == 13 P_T:", np.allclose(G@P_T, 13*P_T), " Mc Mc* P_T:", np.allclose(Mcn@Mcn.conj().T@P_T,13*P_T))
    SMn = S @ Mcn.conj().T          # r Mc*
    Sk = [S[:,perm16[k]] for k in range(7)]
    n=len(S); ch=2000
    n0=0; parts=[]
    for i0 in range(0,n,ch):        # chunks of rows: the pairs passing k = 0, then filtered by k = 1..6 (little memory)
        V = -(SMn[i0:i0+ch] @ Sk[0].T)/16
        ii,jj = np.nonzero(quatok(V))
        n0 += len(ii)
        ii = ii + i0
        for k in range(1,7):
            if len(ii)==0: break
            v = -np.einsum('ij,ij->i', SMn[ii], Sk[k][jj])/16
            m_ = quatok(v)
            ii, jj = ii[m_], jj[m_]
        parts.append(np.stack([ii,jj],axis=1))
    say(name,"pairs passing k=0:",n0)
    good[name] = np.concatenate(parts).astype(np.int64).reshape(-1,2)
    say(name,"seed pairs good at all k:",len(good[name]))
    np.save(os.path.join(out,f"good_{name}.npy"),good[name])
good2 = good["p2"]
SM = S@Mc.conj().T
def vvals(i,j):
    return np.array([-(SM[i]@S[j][perm16[k]])/16 for k in range(7)])

# ---- stage orth
zc=Counter(); Z0={}
for (i,j) in good2:
    v=vvals(i,j); nz=int(np.sum(np.abs(v)<1e-6)); zc[nz]+=1; Z0[(int(i),int(j))]=nz
say("zero counts per (r,c):",dict(zc))
rows=sorted(set(good2[:,0])); cols=sorted(set(good2[:,1]))
say("good pairs",len(good2),"distinct r",len(rows),"distinct c",len(cols))
def orth(i,j):
    # <r_i, tau^k r_j> = 0 for k=0..6
    return all(abs(np.sum(S[i]*np.conj(S[j][perm16[k]])))<1e-9 for k in range(7))
R={}
for i in rows:
    R[i]=[j for j in rows if j>i and orth(i,j)]
say("row seeds",len(rows),"orthogonal row pairs",sum(len(v) for v in R.values()))
Cc={}
for i in cols:
    Cc[i]=[j for j in cols if j>i and orth(i,j)]
say("col seeds",len(cols),"orthogonal col pairs",sum(len(v) for v in Cc.values()))

# ---- stage borders
Cset=set()
for i in cols:
    for j in Cc[i]: Cset.add((int(i),int(j)))
sols=[]
for r1 in rows:
    for r2 in R[r1]:
        cs=[int(c) for c in cols if (int(r1),int(c)) in Z0 and (int(r2),int(c)) in Z0 and Z0[(int(r1),int(c))]+Z0[(int(r2),int(c))]==1]
        for a in range(len(cs)):
            for b in range(a+1,len(cs)):
                c1,c2=cs[a],cs[b]
                if (c1,c2) in Cset or (c2,c1) in Cset:
                    if Z0[(int(r1),c1)]+Z0[(int(r1),c2)]==1 and Z0[(int(r2),c1)]+Z0[(int(r2),c2)]==1:
                        sols.append((int(r1),int(r2),c1,c2))
say("borders (ansatz tau-invariant, seeds):",len(sols))
np.save(os.path.join(out,"borders_p2.npy"),np.array(sols,dtype=np.int64))

# ---- stage check / save: the first border
r1,r2,c1,c2 = sols[0]
def orbit(v): return np.array([v[perm16[k]] for k in range(7)])
Rh = np.concatenate([orbit(S[r1]), orbit(S[r2])],axis=0)       # 14 x 16
Ch = np.concatenate([orbit(S[c1]), orbit(S[c2])],axis=0).T     # 16 x 14
def isq(x): return np.isclose(abs(x),1) and (abs(x.real)<1e-9 or abs(x.imag)<1e-9)
Xb = -Rh@Mc.conj().T@Ch/16
say("R R* = 16 I:", np.allclose(Rh@Rh.conj().T, 16*np.eye(14)), " C* C = 16 I:", np.allclose(Ch.conj().T@Ch, 16*np.eye(14)))
say("halves zero-sum R:", np.allclose(Rh[:,:8].sum(1),0) and np.allclose(Rh[:,8:].sum(1),0), " C:", np.allclose(Ch[:8].sum(0),0) and np.allclose(Ch[8:].sum(0),0))
say("entries quaternary:", all(isq(x) for x in Rh.flatten()), all(isq(x) for x in Ch.flatten()))
say("Xb entries in {0,units}:", all(abs(x)<1e-9 or isq(x) for x in Xb.flatten()))
say("zeros per row:", (np.abs(Xb)<1e-9).sum(1).tolist(), " per col:", (np.abs(Xb)<1e-9).sum(0).tolist())
say("Xb Xb* = 13 I:", np.allclose(Xb@Xb.conj().T, 13*np.eye(14)))
zp = np.abs(Xb)<1e-9
say("zero pattern is a permutation matrix:", bool((zp.sum(0)==1).all() and (zp.sum(1)==1).all()))
def expo(Zm):
    # exponent of i
    return np.rint(np.angle(Zm)/(np.pi/2)).astype(int)%4
js={"Q":"Z2^3 coordinates a=(a0,a1,a2) index 4a0+2a1+a2; tau companion matrix [[0,0,1],[1,0,1],[0,1,0]] acting on column vectors mod 2",
     "compressions(pattern (21,21,6,6))":{"yA":[0,3],"yB":[0,3],"yC":[-1,1],"yD":[-1,1],"meaning":"(y(0), y(a) for a!=0)"},
     "Rhat_exponents_of_i(14x16)":expo(Rh).tolist(),"Chat_exponents_of_i(16x14)":expo(Ch).tolist(),
     "Xb_exponent_or_-1_for_zero(14x14)":[[(-1 if abs(x)<1e-9 else int(expo(np.array([x]))[0])) for x in row] for row in Xb],
     "number_of_tau_invariant_seed_borders":len(sols)}
path = os.path.join(out,"border_k8_pattern2.json")
with open(path,"w") as f: json.dump(js,f)
say("saved",path,"; seeds r1,r2,c1,c2 =",r1,r2,c1,c2)

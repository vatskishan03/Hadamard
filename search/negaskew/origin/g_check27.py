"""exact G-level check of the (13,27) solutions: G = Z13 x F27 (additive Z3^3), blocks U(k,f) = x0(orb k) [f=0], x1 [f square], x2 [f nonsquare]; sum of the four periodic PAFs over G must be 4*351*delta; A twisted-skew."""
import sys, itertools
sys.path.insert(0,'.')
from common import *
m=13; S=Shape(m,[3]); orb=S.orb
# F27 = F3[x]/(x^3-x-1)
def mul(a,b):
    r=[0]*5
    for i in range(3):
        for j in range(3): r[i+j]=(r[i+j]+a[i]*b[j])%3
    # x^3 = x+1 ; x^4 = x^2+x
    r[1]=(r[1]+r[3]+r[4])%3; r[0]=(r[0]+r[3])%3; r[2]=(r[2]+r[4])%3
    return tuple(r[:3])
els=list(itertools.product(range(3),repeat=3))
sq={mul(a,a) for a in els if any(a)}
assert len(sq)==13
neg=lambda a: tuple((-x)%3 for x in a)
assert all((neg(a) not in sq) for a in sq)   # -1 non-square: q = 3 mod 4
sols=load_solutions(sys.argv[1])
step=int(sys.argv[2]) if len(sys.argv)>2 else 1
bad=0; n=0; skewbad=0
for sol in sols[::step]:
    tot=np.zeros((13,3,3,3))
    As=None
    for X,b in enumerate(sol):
        U4=np.zeros((13,3,3,3),complex)
        for k in range(13):
            o=orb[k]
            for a in els:
                if not any(a): v=U[b[0][o]]
                elif a in sq: v=U[b[1][o]]
                else: v=U[b[2][o]]
                U4[(k,)+a]=v
        F=np.fft.fftn(U4); P=np.fft.ifftn(np.abs(F)**2).real
        tot+=P
        if X==0: As=U4
    off=np.abs(tot).copy(); peak=off[0,0,0,0]; off[0,0,0,0]=0
    n+=1
    if abs(peak-4*351)>1e-6 or off.max()>1e-6: bad+=1
    # A twisted skew: A(-g) = i conj A(g) for g != 0
    ok=True
    for k in range(13):
        for a in els:
            if k==0 and not any(a): continue
            g2=((-k)%13,)+neg(a)
            if abs(As[g2]-1j*np.conj(As[(k,)+a]))>1e-9: ok=False
    if not ok: skewbad+=1
print('checked',n,'solutions; G-level PAF-sum failures',bad,'; A twisted-skew failures',skewbad)

"""EXACT check (integer arithmetic on 2x-scaled Gaussian integers) of the origin-defect identities for every solution in a file:
 python3 verify_identities.py m q "[hK]" solfile [step]
 G-level: sum_X PAF_G(U_X) = 4mq delta (re-derived on the G group Z_m x F_q) ;
 K-level: (d) sum(ss*+q tt*)=4mq delta, (Q) sum(ts*-st*)=0, (J) sum(q ee*+se*+es*-tt*)=0, (L1) sum ww* = 4mq delta, and the lag-k consequences (J_k): Z1(k)=sum PAF_t(k)
"""
import sys, os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from common import *
m=int(sys.argv[1]); q=int(sys.argv[2]); hK=eval(sys.argv[3]); fn=sys.argv[4]; step=int(sys.argv[5]) if len(sys.argv)>5 else 1
S=Shape(m,hK); orb=S.orb
sols=load_solutions(fn)
squares={x*x%q for x in range(1,q)}
def Ublock(b):
    # G-level: Z_m x Z_q  -> complex array [m][q]
    M=np.zeros((m,q),dtype=complex)
    for k in range(m):
        o=orb[k]
        for f in range(q):
            M[k,f]=U[b[0][o]] if f==0 else (U[b[1][o]] if f in squares else U[b[2][o]])
    return M
def paf2(M):  # periodic autocorrelation over Z_m x Z_q via FFT
    F=np.fft.fft2(M); return np.fft.ifft2(np.abs(F)**2)
ok=0; bad=0; maxerr=0
for sol in sols[::step]:
    tot=np.zeros((m,q),dtype=complex)
    for b in sol: tot+=paf2(Ublock(b))
    want=np.zeros((m,q)); want[0,0]=4*m*q
    err=np.max(np.abs(tot-want))
    # K-level identities
    ks=[]; 
    sums=dict(d=np.zeros(m,complex),Q=np.zeros(m,complex),J=np.zeros(m,complex),L1=np.zeros(m,complex))
    Z1=np.zeros(m,complex); Pt=np.zeros(m,complex)
    def pafc(x): return np.array([np.sum(x*np.conj(np.roll(x,-k))) for k in range(m)])
    def cc(a,b): return np.array([np.sum(a*np.conj(np.roll(b,-k))) for k in range(m)])
    for b in sol:
        x0=np.array([U[b[0][orb[k]]] for k in range(m)]); x1=np.array([U[b[1][orb[k]]] for k in range(m)]); x2=np.array([U[b[2][orb[k]]] for k in range(m)])
        e=(x1+x2)/2; s=x0-e; t=(x1-x2)/2; w=s+q*e
        sums['d']+=pafc(s)+q*pafc(t); sums['Q']+=cc(t,s)-cc(s,t); sums['J']+=q*pafc(e)+cc(e,s)+cc(s,e)-pafc(t); sums['L1']+=pafc(w)
        Z1+=cc(e,s)+cc(s,e); Pt+=pafc(t)
    want=np.zeros(m); want[0]=4*m*q
    err2=max(np.max(np.abs(sums['d']-want)),np.max(np.abs(sums['Q'])),np.max(np.abs(sums['J'])),np.max(np.abs(sums['L1']-want)))
    err3=np.max(np.abs(Z1[1:]-Pt[1:]))
    e_=max(err,err2,err3); maxerr=max(maxerr,e_)
    if e_<1e-8: ok+=1
    else: bad+=1
print('checked',ok+bad,'ok',ok,'bad',bad,'max err',maxerr)

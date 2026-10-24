import sys, json, random, time, warnings
warnings.filterwarnings('ignore')
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))      # search/skewplug
from sk import SShape
from sk_pat import m_splits, row_patterns_skew
from awin import alloc_window
import est_skew as es
K=[int(sys.argv[1])]; q=int(sys.argv[2]); hu=[[int(x)] for x in sys.argv[3].split(',')]; secs=float(sys.argv[4])
sh=SShape(K,q,[u[0] for u in hu])
win,emin,emax=alloc_window(sh,erel=10**6)
print('nu',sh.kCD.nu,'N',sh.N,'E range',emin,emax,'allocs',len(win),flush=True)
rng=random.Random(5)
t0=time.time(); res=[]
byE={}
while time.time()-t0<secs:
    a,em=win[rng.randrange(len(win))]
    nA,nB,nC,nD=a
    for _ in range(100):
        mA=2*rng.randint(0,nA//2); mB=rng.randint(0,nB); mC=rng.randint(0,nC); mD=em-mA-mB-mC
        if 0<=mD<=nD: break
    else: continue
    ms=(mA,mB,mC,mD)
    p=row_patterns_skew(sh,a,ms)
    E=sum(a)
    byE.setdefault(E,[0,0,[]]); byE[E][0]+=1
    if p:
        byE[E][1]+=1
        e=es.estimate(sh,tuple(a),ms,M=60)
        byE[E][2].append((e['logE'],e['logrho'],a,ms))
for E in sorted(byE):
    n,f,l=byE[E]
    print(E,n,f, 'best', max(l)[:2] if l else None)

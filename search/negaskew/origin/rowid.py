import sys
sys.path.insert(0,'.')
from common import *
from collections import Counter
m=int(sys.argv[1]); q=int(sys.argv[2]); hK=eval(sys.argv[3]); fn=sys.argv[4]
S=Shape(m,hK); orb=np.array(S.orb); ell=S.ell; h=ell[1]; p=(S.nu-1)//2
sols=load_solutions(fn); Uarr=np.array(U)
c=Counter()
for sol in sols[::max(1,len(sols)//3000)]:
    ds=0; dt=0; hs=0; ht=0
    for b in sol:
        x0=Uarr[np.array(b[0])[orb]]; x1=Uarr[np.array(b[1])[orb]]; x2=Uarr[np.array(b[2])[orb]]
        e=(x1+x2)/2; s=x0-e; t=(x1-x2)/2
        Rs=s.sum(); Rt=t.sum(); Ss=s[0]; Tt=t[0]
        ds+=abs(Rs)**2-abs(Ss)**2; dt+=abs(Rt)**2-abs(Tt)**2
    c[(round(ds,6),round(dt,6))]+=1
print('m',m,'h',h,'p',p,'4(m-1)=',4*(m-1),'(sum(|R^s|^2-|S|^2), sum(|R^t|^2-|tau|^2)) counts:',dict(c))

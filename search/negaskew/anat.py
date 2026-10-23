"""Anatomy of a plug-in solution of qp: for each block the type of every K-orbit (E: x2 = -x1, no defect; F: full defect x2 = x1;
M+, M-: half defects x2 = +-i x1) with its defect cost (q-1)|e|^2 + 2 Re(x0 conj e), the weight W of each block and the total, which is 4|K|.
usage: python3 anat.py m q "[H_K generators]" "SOL line (blocks x0|x1|x2)"
Both hits have every defect at the origin k = 0 of K: v = 299, K = Z_13, q = 23, H_K = <3>:
  python3 anat.py 13 23 "[3]" "03120|31320|23102 00313|10322|32100 01233|20130|22312 01022|30132|32310"
gives W = 10 (A, half defect), 0 (B), 20 (C), 22 (D) and sum W = 52 = 4|K|; v = 589 (K = Z_19, q = 31, H_K = <7>): A 16, B 32, D 28 (NOTES.md section 89)."""
import sys, json
from qshape import kshape
def anal(m,q,hK,sol):
    sh=kshape([m],hK); ell=sh['ell']
    U=[1,1j,-1,-1j]
    tot=0
    for X,b in zip('ABCD',sol.split()):
        x0,x1,x2=b.split('|'); line=[]; W=0
        for o in range(sh['nu']):
            a0,a1,a2=int(x0[o]),int(x1[o]),int(x2[o])
            d=(a2-a1)%4
            typ={2:'E',0:'F',1:'M+',3:'M-'}[d]
            e=(U[a1]+U[a2])/2; x=U[a0]
            c=(q-1)*abs(e)**2+2*(x*e.conjugate()).real
            W+=ell[o]*round(c)
            line.append(f"{typ}{round(c) if typ!='E' else ''}")
        tot+=W
        print(X,'W=',W,' '.join(line))
    print('sum W',tot,'4|K|',4*m,'ells',ell)
anal(int(sys.argv[1]),int(sys.argv[2]),json.loads(sys.argv[3]),sys.argv[4])

"""H(8(g1*m + n)) from a DCS (a, b) of length m free of +-is (Kamali-Kharaghani, Lemma 3, Theorem 7, Lemmas 2 and 1).
usage: python3 kk_big.py A B G1 N [OUT.npy]   (A, B comma-separated tokens 1, -1, i, -i, s, -s; G1, N Golay lengths)
Checks H H^T = nI exactly (float32 products of +-1 matrices are exact below 2^24), prints the SHA-256 of the CSV
(rows of 1 and -1, comma-separated) and saves H to OUT.npy if given."""
import numpy as np, sys, os, hashlib, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kk_build import el, golay, dcs_times_golay, I2, S
a=sys.argv[1].split(','); b=sys.argv[2].split(','); g1=int(sys.argv[3]); n=int(sys.argv[4])
out=sys.argv[5] if len(sys.argv)>5 else None
t0=time.time()
a,b=dcs_times_golay(a,b,g1)
c,d=golay(n); m=len(a); L=m+n
def rowseq(x,y,sgn):
    return [el(t) for t in x]+[sgn*v*I2 for v in y]
first={'A':rowseq(a,c,1),'B':rowseq(a,c,-1),'C':rowseq(b,d,1),'D':rowseq(b,d,-1)}
# H = rho(M)(I (x) F); build directly in int8, block (r,c) of the 4x4 array
F=np.array([[1,1],[1,-1]],dtype=np.int8)
N=8*L
H=np.zeros((N,N),dtype=np.int8)
idx=np.arange(L)
def circ_block(seq,transpose=False,times_Rs=False,sign=1):
    # returns 2L x 2L int8 of sign * X (or X^T) (times Rs)
    X=np.zeros((2*L,2*L),dtype=np.int8)
    for k in range(L):
        blk=seq[k]
        # X[r, c] = seq[(c-r)%L]
        r=idx; cc=(idx+k)%L
        if transpose: r,cc=cc,r; blk=blk.T
        for i in range(2):
            for j in range(2):
                X[2*r+i,2*cc+j]=blk[i,j]
    if times_Rs:
        Y=np.zeros_like(X)
        # (X Rs)[:, block L-1-c] = X[:, block c] @ S
        for cb in range(L):
            Y[:,2*(L-1-cb):2*(L-1-cb)+2]=X[:,2*cb:2*cb+2]@S.astype(np.int8)
        X=Y
    return sign*X
spec=[[('A',0,0,1),('B',0,1,-1),('C',0,1,-1),('D',0,1,-1)],
      [('B',0,1,1),('A',0,0,1),('D',1,1,-1),('C',1,1,1)],
      [('C',0,1,1),('D',1,1,1),('A',0,0,1),('B',1,1,-1)],
      [('D',0,1,1),('C',1,1,-1),('B',1,1,1),('A',0,0,1)]]
for R in range(4):
    for Cc in range(4):
        nm,tr,rs,sg=spec[R][Cc]
        X=circ_block(first[nm],bool(tr),bool(rs),sg)
        # times (I (x) F)
        Xf=np.empty_like(X)
        Xf[:,0::2]=X[:,0::2]+X[:,1::2]; Xf[:,1::2]=X[:,0::2]-X[:,1::2]
        H[2*L*R:2*L*(R+1),2*L*Cc:2*L*(Cc+1)]=Xf
print('built',N,'in',round(time.time()-t0,1),'s; entries +-1:',bool(np.all(np.abs(H)==1)))
Hf=H.astype(np.float32)
G=Hf@Hf.T
del Hf
np.fill_diagonal(G, G.diagonal() - N)
ok = not np.any(G)
del G
print('H H^T = %d I (float32 BLAS, exact for |entries|<2^24):'%N, ok, 'time',round(time.time()-t0,1))
s=hashlib.sha256()
for row in H:
    s.update((','.join('1' if v>0 else '-1' for v in row)+'\n').encode())
print('sha256 csv',s.hexdigest())
if out:
    np.save(out,H)

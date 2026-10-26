"""The quaternary conference matrix K (order 478) of the n = 955 core from a solution of s1_prod.py and the border.

usage: python3 assemble_K.py SOLUTION.json [--border border_k8_pattern2.json] [--out DIR]
SOLUTION.json: four +-1 sequences A, B, C, D on G = Z2^3 x Z29 (index ((e0*2+e1)*2+e2)*29+z); the border file defaults to
data/border_k8_pattern2.json (border8.py regenerates it). Builds Kc = [[X, Y], [Y*, -X*]] (Szoellosi's array over G with
X = alpha X_A + conj(alpha) X_B, Y = alpha X_C + conj(alpha) X_D R, alpha = (1 + i)/2) and K = [[Xb, R], [R', Kc]] with the
border of complex-hadamard/TBORDERED.md (rows of Rhat and columns of Chat lifted through G -> Q = Z2^3), prints the exact
checks K K* = 477 I, one zero per row and column, entries in {1, i, -1, -i}, and saves OUT/K478.npy (complex128).
"""
import argparse, json, os, sys, numpy as np, itertools
import tbio

ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("solution")
ap.add_argument("--border", default=os.path.join(tbio.DATA, "border_k8_pattern2.json"))
ap.add_argument("--out", default=None, help="output directory (default: current directory)")
args = ap.parse_args()
out = tbio.outdir(args.out)

sol=json.load(open(args.solution))
A,B,C,D=[np.array(sol[k],dtype=np.int64) for k in "ABCD"]
mods=(2,2,2,29); m=232
coords=np.array(list(itertools.product(*[range(x) for x in mods])))
mo=np.array(mods)
def idx(c): return np.ravel_multi_index(tuple((c%mo).T),mods)
diff=idx(coords[None,:,:]-coords[:,None,:])       # diff[g,h] = index of h-g
def dev(X): return X[diff]                         # X_Z[g,h]=Z(h-g)
negp=idx(-coords)
R=np.zeros((m,m)); R[np.arange(m),negp]=1
al=(1+1j)/2
XA,XB,XC,XD=[dev(v).astype(float) for v in (A,B,C,D)]
X=al*XA+np.conj(al)*XB
Y=al*XC+np.conj(al)*(XD@R)
Kc=np.block([[X,Y],[Y.conj().T,-X.conj().T]])
# border
bd=json.load(open(args.border))
Rh=1j**np.array(bd["Rhat_exponents_of_i(14x16)"]); Ch=1j**np.array(bd["Chat_exponents_of_i(16x14)"])
Xb=np.array([[0 if e==-1 else 1j**e for e in row] for row in bd["Xb_exponent_or_-1_for_zero(14x14)"]])
phi=coords[:,:3]%2; qi=phi[:,0]*4+phi[:,1]*2+phi[:,2]      # index in Q, 4a0+2a1+a2
Rm=np.zeros((14,2*m),dtype=complex); Cm=np.zeros((2*m,14),dtype=complex)
for i in (0,1):
    Rm[:,i*m:(i+1)*m]=Rh[:,i*8+qi]
    Cm[i*m:(i+1)*m,:]=Ch[i*8+qi,:]
K=np.block([[Xb,Rm],[Cm,Kc]])
n=K.shape[0]
G=K@K.conj().T
print("order",n,"K K* == (n-1) I :",np.allclose(G,(n-1)*np.eye(n)))
zeros=np.abs(K)<1e-9
print("zeros per row all 1:",(zeros.sum(1)==1).all()," per col all 1:",(zeros.sum(0)==1).all())
nz=K[~zeros]
print("nonzero entries units in {1,i,-1,-i}:",bool(np.all(np.isclose(np.abs(nz),1))) and bool(np.all((np.abs(nz.real)<1e-9)|(np.abs(nz.imag)<1e-9))))
np.save(os.path.join(out,'K478.npy'),K)
print("wrote", os.path.join(out,'K478.npy'))

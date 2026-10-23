"""Pairing test of octic quadruples: is a solution (A, B, C, D) made of two theta-pairs (X, m theta(X)) (theta: x -> -x, m a monomial +-x^s, with adjoints),
i.e. the image of an order-8v quadruple with a centre-palindromic block?  usage: python3 dbl.py SOLFILE n [limit]   (SOL lines of exsmall: four integers = bit
masks of the sequences of length n = 4v; e.g. ./exsmall 12 1000000 0 > runs/ex12.out; python3 dbl.py runs/ex12.out 12 3000 -> tested 3000, doubled-type 168 (the first 3000 solutions; all 148,848: 10,368 = 7.0%, ./exsmall 12 0 1)).
module: mono_images, theta, adj, paired, pair_any, is_doubled"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, sys
from mkcert import bits_to_seq
def mono_images(Y):
    n=len(Y); out=[]
    cur=np.array(Y)
    for s in range(n):
        out.append(cur.copy()); out.append(-cur)
        nxt=np.empty_like(cur); nxt[1:]=cur[:-1]; nxt[0]=-cur[-1]; cur=nxt
    return out
def theta(Y):  # x -> -x
    n=len(Y); return np.array([Y[k]*(-1)**k for k in range(n)])
def adj(Y):
    n=len(Y); o=np.empty_like(Y); o[0]=Y[0]; o[1:]=-Y[:0:-1]; return o
def paired(X,Y):
    """is Y = monomial * theta(X) or monomial*theta(X)^* ?  (X given as the first of the pair)"""
    tX=theta(X); cands=[tX,adj(tX)]
    for c in cands:
        for im in mono_images(c):
            if np.array_equal(im,Y): return True
    return False
def pair_any(X,Y):
    # symmetric test with the allowed variants (adjoint of Y and X)
    return paired(X,Y) or paired(X,adj(Y)) or paired(adj(X),Y)
def is_doubled(S):
    A,B,C,D=S
    for i,(X,Y,Z) in enumerate(((B,C,D),(C,B,D),(D,B,C))):
        if pair_any(A,X) and pair_any(Y,Z): return True
    return False
if __name__=="__main__":
    fn=sys.argv[1]; n=int(sys.argv[2]); lim=int(sys.argv[3]) if len(sys.argv)>3 else 10**9
    tot=0; dbl=0; ex=[]
    for l in open(fn):
        if not l.startswith('SOL'): continue
        S=[np.array(bits_to_seq(int(x),n)) for x in l.split()[1:]]
        tot+=1
        if is_doubled(S): dbl+=1
        elif len(ex)<3: ex.append(S)
        if tot>=lim: break
    print('solutions tested',tot,' doubled-type (theta-paired up to monomials/adjoints):',dbl,' not doubled:',tot-dbl)

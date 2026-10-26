# Independent DCS checker in the 2x2 real-matrix model of the dihedral signed group.
# 1->I, i->J=[[0,-1],[1,0]], s->S=[[1,0],[0,-1]], is = J@S. Involution = transpose.
# Checks (i) f_a,f_b,f_a*,f_b* commute as matrix Laurent polynomials, (ii) F_a F_a^T(1/t) + F_b F_b^T(1/t) = 2n I.
import numpy as np, sys, json
I=np.eye(2,dtype=int); J=np.array([[0,-1],[1,0]]); S=np.array([[1,0],[0,-1]])
EL={'1':I,'i':J,'s':S,'is':J@S}
def parse(tok):
    sg=-1 if tok.startswith('-') else 1; return sg*EL[tok.lstrip('-')]
def poly(seq):  # list of 2x2 coefficient matrices, degree index j
    return np.array([parse(t) for t in seq])
def star(P):    # involution: transpose coefficients, t -> 1/t  (store as dict deg->mat)
    return {-j:P[j].T for j in range(len(P))}
def asdict(P): return {j:P[j] for j in range(len(P))}
def mul(X,Y):
    Z={}
    for i,A in X.items():
        for j,B in Y.items():
            Z[i+j]=Z.get(i+j,np.zeros((2,2),int))+A@B
    return Z
def eq(X,Y):
    ks=set(X)|set(Y)
    return all(np.array_equal(X.get(k,np.zeros((2,2),int)),Y.get(k,np.zeros((2,2),int))) for k in ks)
def check(a,b):
    n=len(a); assert len(b)==n
    A,B=poly(a),poly(b)
    fs=[asdict(A),asdict(B),star(A),star(B)]
    comm=all(eq(mul(fs[x],fs[y]),mul(fs[y],fs[x])) for x in range(4) for y in range(x+1,4))
    S2=mul(asdict(A),star(A)); T=mul(asdict(B),star(B))
    tot={k:S2.get(k,0)+T.get(k,0) for k in set(S2)|set(T)}
    golay=all(np.array_equal(v,(2*n*I if k==0 else 0*I)) for k,v in tot.items())
    noIS=all(t.lstrip('-')!='is' for t in a+b)
    return comm,golay,noIS
if __name__=='__main__':
    a=sys.argv[1].split(','); b=sys.argv[2].split(',')
    print(check(a,b))

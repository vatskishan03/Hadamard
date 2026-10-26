# Kamali-Kharaghani chain: DCS (a,b) length m (free of is) + Golay (c,d) length n ->
# D-circulants A=circ(a,c), B=circ(a,-c), C=circ(b,d), D=circ(b,-d) of order L=m+n ->
# Craigen-Goethals-Seidel array (dihedral, order 4L) -> real H(8L) via rho(M)(I (x) F).
import numpy as np, sys
I2=np.eye(2,dtype=np.int64); J=np.array([[0,-1],[1,0]]); S=np.array([[1,0],[0,-1]])
EL={'1':I2,'i':J,'s':S,'is':J@S}
def el(tok):
    sg=-1 if tok.startswith('-') else 1; return sg*EL[tok.lstrip('-')]
def golay(n):
    """A real Golay pair of length n = 2^a 10^b 26^c, built from the base pairs of lengths 2, 10 and 26 by the
    Kronecker product rule e = (c+d)/2 (x) a + (c-d)/2 (x) b, f = rev((d-c)/2) (x) a + rev((c+d)/2) (x) b."""
    base={2:([1,1],[1,-1]),10:([1,1,-1,1,-1,1,-1,-1,1,1],[1,1,-1,1,1,1,1,1,-1,-1]),
          26:([1,1,1,1,-1,1,1,-1,-1,1,-1,1,-1,1,-1,-1,1,-1,1,1,1,-1,-1,1,1,1],
              [1,1,1,1,-1,1,1,-1,-1,1,-1,1,1,1,1,1,-1,1,-1,-1,-1,1,1,-1,-1,-1])}
    def ok(a,b):
        m=len(a)
        return all(int(a[:m-k]@a[k:]+b[:m-k]@b[k:])==(2*m if k==0 else 0) for k in range(m))
    a=b=np.array([1]); x=n
    for g in (26,10,2):
        while x%g==0:
            c,d=map(np.array,base[g])
            a,b=np.kron((c+d)//2,a)+np.kron((c-d)//2,b),np.kron(((d-c)//2)[::-1],a)+np.kron(((c+d)//2)[::-1],b)
            x//=g
    if x!=1 or not ok(a,b):
        raise ValueError(f"no Golay pair of length {n} from this rule")
    return [int(v) for v in a],[int(v) for v in b]
def dcirc(seq):  # D-circulant of first row seq (list of 2x2 matrices), as 2L x 2L real block matrix
    L=len(seq); M=np.zeros((2*L,2*L),dtype=np.int64)
    for r in range(L):
        for c in range(L):
            M[2*r:2*r+2,2*c:2*c+2]=seq[(c-r)%L]
    return M
def blockT(M):  # involution on block matrices: transpose blocks and block positions = plain transpose
    return M.T
def build(a,b,c,d,variant=0):
    m=len(a); n=len(c); L=m+n
    A=[el(t) for t in a]+[ci*I2 for ci in c]
    B=[el(t) for t in a]+[-ci*I2 for ci in c]
    C=[el(t) for t in b]+[di*I2 for di in d]
    D=[el(t) for t in b]+[-di*I2 for di in d]
    A,B,C,D=map(dcirc,(A,B,C,D))
    R=np.zeros((2*L,2*L),dtype=np.int64)
    for r in range(L): R[2*r:2*r+2,2*(L-1-r):2*(L-1-r)+2]=S   # R s : back diagonal times s
    Rs=R
    T=lambda X: X.T
    rows=[[A,-B@Rs,-C@Rs,-D@Rs],[B@Rs,A,-T(D)@Rs,T(C)@Rs],[C@Rs,T(D)@Rs,A,-T(B)@Rs],[D@Rs,-T(C)@Rs,T(B)@Rs,A]]
    M=np.block(rows)
    return M
if __name__=='__main__':
    a=sys.argv[1].split(','); b=sys.argv[2].split(','); n=int(sys.argv[3])
    c,d=golay(n)
    M=build(a,b,c,d)
    N=M.shape[0]
    G=M@M.T
    ok=np.array_equal(G,(N//2)*np.eye(N,dtype=np.int64))
    print('dihedral order',N//2,'rho(M)rho(M)^T = (N/2) I:',ok)
    F=np.array([[1,1],[1,-1]])
    H=M@np.kron(np.eye(N//2,dtype=np.int64),F)
    print('H order',H.shape[0],'entries +-1:',bool(np.all(np.abs(H)==1)),'HH^T=nI:',np.array_equal(H@H.T,H.shape[0]*np.eye(H.shape[0],dtype=np.int64)))

def neg(tok):
    return tok[1:] if tok.startswith('-') else '-' + tok

def dcs_times_golay(a, b, g1):
    """Lemma 3 of Kamali-Kharaghani: a DCS (a, b) of length m (token lists) and a Golay pair (c, d) of length g1 give
    the DCS e = (c+d)/2 (x) a + (c-d)/2 (x) b, h = rev((c-d)/2) (x) a - rev((c+d)/2) (x) b of length m g1."""
    c, d = golay(g1)
    p = [(x + y) // 2 for x, y in zip(c, d)]
    q = [(x - y) // 2 for x, y in zip(c, d)]
    e, h = [], []
    for k in range(g1):
        for r in range(len(a)):
            e.append(a[r] if p[k] == 1 else neg(a[r]) if p[k] == -1 else b[r] if q[k] == 1 else neg(b[r]))
            pk, qk = p[g1 - 1 - k], q[g1 - 1 - k]
            h.append(a[r] if qk == 1 else neg(a[r]) if qk == -1 else neg(b[r]) if pk == 1 else b[r])
    return e, h

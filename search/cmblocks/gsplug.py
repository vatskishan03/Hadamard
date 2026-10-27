import numpy as np, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import cmb

def real_blocks(N):
    """B_rc = [[X, Y], [-Y, X]], X = Hr + Hi, Y = Hr - Hi (order 2q, int8), from complex quaternary blocks"""
    out = [[None] * 4 for _ in range(4)]
    for r in range(4):
        for c in range(4):
            Z = N[r][c]
            Hr, Hi = np.rint(Z.real).astype(np.int8), np.rint(Z.imag).astype(np.int8)
            X, Y = Hr + Hi, Hr - Hi
            out[r][c] = np.block([[X, Y], [-Y, X]]).astype(np.int8)
    return out

def Rprime(q):
    """signed permutation R' = diag(1,-1) (x) R_B on the order-2q index (eta, l): returns (perm, sign) with R'[i, perm[i]] = sign[i]"""
    Rb = cmb.reversal_perm(q)
    perm = np.concatenate([Rb, q + Rb]); sign = np.concatenate([np.ones(q, np.int8), -np.ones(q, np.int8)])
    return perm, sign

def sperm_dense(perm, sign):
    n = len(perm); M = np.zeros((n, n), np.int64); M[np.arange(n), perm] = sign
    return M

def check_real_blocks(Bk, q):
    perm, sign = Rprime(q)
    Rm = sperm_dense(perm, sign)
    flat = [Bk[r][c].astype(np.int64) for r in range(4) for c in range(4)]
    ok_rt = all(np.array_equal(Rm @ x @ Rm, x.T) for x in flat)
    ok_comm = all(np.array_equal(x @ y, y @ x) for i, x in enumerate(flat) for y in flat[i:])
    n = 2 * q
    ok_orth = all(np.array_equal(sum(Bk[r][c].astype(np.int64) @ Bk[s][c].astype(np.int64).T for c in range(4)),
                                 (8 * q if r == s else 0) * np.eye(n, dtype=np.int64)) for r in range(4) for s in range(4))
    return ok_rt, ok_comm, ok_orth

def plug_A(Bk, T):
    """A_j = sum_k T_k (x) B_kj, order t * 2q (int8)"""
    Ct = [cmb.circ_int(r).astype(np.int8) for r in T]
    return [sum(np.kron(Ct[k], Bk[k][j]) for k in range(4)).astype(np.int8) for j in range(4)]

def gs_R(q, t):
    """R'' = R_t (x) R' as (perm, sign) on index (s, eta, l)"""
    perm, sign = Rprime(q)
    n2 = 2 * q
    P = np.empty(t * n2, np.int64); S = np.empty(t * n2, np.int8)
    for s in range(t):
        P[s * n2:(s + 1) * n2] = ((-s) % t) * n2 + perm
        S[s * n2:(s + 1) * n2] = sign
    return P, S

def mulR(X, P, S):
    """X @ R'' for R''[i, P[i]] = S[i]:  (X R'')[:, P[i]] = X[:, i] S[i]"""
    Y = np.empty_like(X)
    Y[:, P] = X * S[None, :]
    return Y

def gs_array(A, P, S):
    n = A[0].shape[0]
    AR = [mulR(x, P, S) for x in A]
    ATR = [mulR(np.ascontiguousarray(x.T), P, S) for x in A]
    rows = [[A[0], AR[1], AR[2], AR[3]], [-AR[1], A[0], ATR[3], -ATR[2]],
            [-AR[2], -ATR[3], A[0], ATR[1]], [-AR[3], ATR[2], -ATR[1], A[0]]]
    G = np.empty((4 * n, 4 * n), np.int8)
    for i in range(4):
        for j in range(4):
            G[i * n:(i + 1) * n, j * n:(j + 1) * n] = rows[i][j]
    return G

def is_hadamard(G):
    F = G.astype(np.float64)
    return np.array_equal(F @ F.T, G.shape[0] * np.eye(G.shape[0]))

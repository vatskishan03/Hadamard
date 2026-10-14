"""Bordered T-matrices -> OD(4t; t,t,t,t) (Baumert-Hall array) -> H(4tw) with Williamson matrices.

Bordered T-matrices of order t: T_1..T_4 : G -> {0,+-1}, G abelian of order v = t - 1, supports partition G,
sum_k P_{T_k}(s) = -1 for s != 0, row sums (+-1, 0, 0, 0).
Core: Cooper-Wallis blocks A..D (linear in the variables a,b,c,d) in the Goethals-Seidel array over G.
Border: H = [[-Q, Q (x) 1^T], [GS(c) (x) 1, GS(A,B,C,D)]], Q = quaternion OD(4;1,1,1,1), c = row sums of A..D."""
import os
import itertools, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import spence_lib as L

def gmat(G, f):                     # M[g, h] = f(h - g)
    idx = G.index(G.coords[None, :, :] - G.coords[:, None, :])
    return f[idx]

def blocks(T, x):
    a, b, c, d = x; T1, T2, T3, T4 = T
    A = a*T1 + b*T2 + c*T3 + d*T4
    B = -b*T1 + a*T2 + d*T3 - c*T4
    C = -c*T1 - d*T2 + a*T3 + b*T4
    D = -d*T1 + c*T2 - b*T3 + a*T4
    return A, B, C, D

def od_matrix(G, T, x):
    A, B, C, D = [gmat(G, s) for s in blocks(T, x)]
    v = G.v
    R = np.zeros((v, v)); R[np.arange(v), G.neg()] = 1
    core = np.block([[A, B @ R, C @ R, D @ R],
                     [-B @ R, A, D.T @ R, -C.T @ R],
                     [-C @ R, -D.T @ R, A, B.T @ R],
                     [-D @ R, C.T @ R, -B.T @ R, A]])
    cA, cB, cC, cD = [s.sum() for s in blocks(T, x)]
    gs = np.array([[cA, cB, cC, cD], [-cB, cA, cD, -cC], [-cC, -cD, cA, cB], [-cD, cC, -cB, cA]])
    a, b, c, d = x
    Q = np.array([[a, b, c, d], [-b, a, -d, c], [-c, d, a, -b], [-d, -c, b, a]])
    one = np.ones((1, v))
    return np.block([[-Q, np.kron(Q, one)], [np.kron(gs, one.T), core]])

def is_bordered_T(G, T):
    T = np.array(T)
    if not np.all(np.abs(T).sum(0) == 1): return False
    P = sum(np.array([t @ t[G.translate(g)] for g in range(G.v)]) for t in T)
    rs = sorted(abs(int(t.sum())) for t in T)
    return P[0] == G.v and np.all(P[1:] == -1) and rs == [0, 0, 0, 1]

def find_small(mods):
    G = L.AbGroup(mods); v = G.v; sols = []
    for cfg in itertools.product(range(8), repeat=v):
        T = np.zeros((4, v))
        for g, s in enumerate(cfg): T[s % 4, g] = 1 if s < 4 else -1
        if is_bordered_T(G, T): sols.append(T)
    return G, sols

if __name__ == '__main__':
    rng = np.random.default_rng(0)
    for mods in ([4], [2, 2], [8], [2, 4], [2, 2, 2]):
        G, sols = find_small(mods)
        print(f"G=Z{mods} (t={G.v+1}): {len(sols)} bordered T-matrix configurations", flush=True)
        if not sols: continue
        T = sols[0]
        ok = True
        for _ in range(5):
            x = rng.standard_normal(4)
            H = od_matrix(G, T, x)
            ok &= np.allclose(H @ H.T, (G.v + 1) * (x @ x) * np.eye(H.shape[0]))
        Hpm = od_matrix(G, T, np.ones(4))
        print(f"   OD(4t; t^4) identity holds for random real x: {ok}; entries at x=1 are +-1: {set(np.unique(Hpm))}", flush=True)

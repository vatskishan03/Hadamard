"""573 Spence-spin structure: G = Z_2^2 x Z_143, H = <27> (order 5: QR on Z_11, trivial on Z_13, trivial on Z_2^2),
mu = (alpha, u), alpha of order 3 on Z_2^2, u = 3 or 9 on Z_143 (mu^3 = (1,u^3) in H)."""
import numpy as np
from ss_core import Grp, Spin
ALPHA = np.array([[0, 1], [1, 1]])
def aut_matrix(sig, u):
    M = np.zeros((3, 3), dtype=np.int64); M[:2, :2] = sig; M[2, 2] = u
    return M
def build(u=3, m=143, hgen=27, alpha_pow=1):
    G = Grp([2, 2, m])
    sig = np.eye(2, dtype=np.int64)
    for _ in range(alpha_pow): sig = (sig @ ALPHA) % 2
    Hp = [G.aut(aut_matrix(np.eye(2, dtype=np.int64), hgen))]
    mu = G.aut(aut_matrix(sig, u))
    return Spin(G, Hp, mu)
if __name__ == '__main__':
    for u in (3, 9):
        S = build(u)
        print('u=%d D=%d kK=%d kH=%d trivial_spin=%s |szK|=%s' % (u, S.D, S.kK, S.kH, S.trivial_spin, sorted(set(S.szK.tolist()))))
        print('  H orbit sizes', dict(zip(*np.unique(S.szH, return_counts=True))), ' K orbit sizes', dict(zip(*np.unique(S.szK, return_counts=True))))
        print('  class sizes', sorted(S.csize.tolist()))
        print('  pool A (row sum 2) 10^%.3f, pool B (row sum 0) 10^%.3f' % (np.log10(S.pool('A', 2)), np.log10(S.pool('B', 0))))
        # H = {1,14,27,53,92}?
        G = S.G
        e = np.flatnonzero((G.coords[:, 0] == 0) & (G.coords[:, 1] == 0))   # copy of Z_143 : index of element m -> coords (0,0,m)
        Hm = sorted({int(G.coords[S.Hperms[0][ [i for i in range(G.v) if tuple(G.coords[i]) == (0,0,1)][0] ]][2])} ) if False else None
        # H elements as multiplier group on Z_143
        x = 1; Hs = []
        for _ in range(5): Hs.append(x); x = x * 27 % 143
        print('  H =', sorted(Hs), ' mu on Z_143 = x%d (order %d)' % (u, next(k for k in range(1, 100) if pow(u, k, 143) == 1)))

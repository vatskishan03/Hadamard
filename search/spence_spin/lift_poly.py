"""Exact polynomial decomposition of the 27 class keys of a lifted (A,B) as function of the Legendre-fibre signs.
X(s) = X0 + D s  (affine, exact);  key_c(sA,sB) = KA_c(sA) + KB_c(sB),
KA_c(s) = cA_c + linA_c . s + sum_{i<j} QA_c[i,j] s_i s_j  (same for B), all integers.
Also: exact linear consequences (integer combinations w of classes whose quadratic parts cancel identically)."""
import numpy as np, sys
from fractions import Fraction
from struct573 import build
from levels import make
from liftlab import Lifter, triple

def base_and_D(lf, a, b):
    """returns LA, LB, X0A, X0B (s=0 lifts), DA (kK x nA), DB (kH x nB)"""
    LB, LA = lf.legendre_points(a, b)
    S = lf.S
    def lift0(sA, sB):
        XB = np.zeros(S.kH, dtype=np.int64); XA = np.zeros(S.kK, dtype=np.int64)
        tb = {q: triple(b[q], sB.get(q, 0)) if abs(b[q]) == 1 else triple(b[q], 1) for q in range(len(b))}
        ta = {R: triple(a[R], sA.get(R, 0)) if abs(a[R]) == 1 else triple(a[R], 1) for R in range(len(a))}
        for o in range(S.kH): XB[o] = tb[int(lf.qB[o])][int(lf.cB[o])]
        for o in range(S.kK): XA[o] = ta[int(lf.qA[o])][int(lf.cA[o])]
        return XA, XB
    # NB triple(+-1, 0) = (+-1, 0, 0) (s=0): affine base
    X0A, X0B = lift0({}, {})
    DA = np.zeros((S.kK, len(LA)), dtype=np.int64); DB = np.zeros((S.kH, len(LB)), dtype=np.int64)
    for i, R in enumerate(LA):
        XA, _ = lift0({R: 1}, {}); DA[:, i] = XA - X0A
    for i, q in enumerate(LB):
        _, XB = lift0({}, {q: 1}); DB[:, i] = XB - X0B
    return LA, LB, X0A, X0B, DA, DB

def poly(M, X0, D):
    """M: (D x k x k) int class matrices.  returns const (D,), lin (D x n), Q (D x n x n) upper-triangular pair coefficient (i<j)"""
    nC = M.shape[0]; n = D.shape[1]
    const = np.zeros(nC, dtype=np.int64); lin = np.zeros((nC, n), dtype=np.int64); Q = np.zeros((nC, n, n), dtype=np.int64)
    for c in range(nC):
        Mc = M[c].astype(np.int64)
        c0 = int(X0 @ Mc @ X0)
        l = D.T @ (Mc + Mc.T) @ X0
        N = D.T @ Mc @ D
        diag = int(np.trace(N))
        Qc = np.triu(N + N.T, 1)
        const[c] = c0 + diag; lin[c] = l; Q[c] = Qc
    return const, lin, Q

class Poly:
    def __init__(self, S, LQ, a, b):
        self.S, self.LQ, self.a, self.b = S, LQ, np.asarray(a, dtype=np.int64), np.asarray(b, dtype=np.int64)
        lf = self.lf = Lifter(S, LQ)
        self.LA, self.LB, self.X0A, self.X0B, self.DA, self.DB = base_and_D(lf, self.a, self.b)
        self.nA, self.nB = len(self.LA), len(self.LB)
        self.cA, self.linA, self.QA = poly(S.MA, self.X0A, self.DA)
        self.cB, self.linB, self.QB = poly(S.MB, self.X0B, self.DB)
    def keys_exact(self, sA, sB):
        """exact class keys via the original (non-polynomial) route: sA,sB in +-1 arrays aligned with LA, LB"""
        XA = self.X0A + self.DA @ np.asarray(sA, dtype=np.int64); XB = self.X0B + self.DB @ np.asarray(sB, dtype=np.int64)
        return np.rint(self.S.valsA(XA[None].astype(float))).astype(np.int64)[0] + np.rint(self.S.valsB(XB[None].astype(float))).astype(np.int64)[0], XA, XB
    def keys_poly(self, sA, sB):
        sA = np.asarray(sA, dtype=np.int64); sB = np.asarray(sB, dtype=np.int64)
        kA = self.cA + self.linA @ sA + np.einsum('i,cij,j->c', sA, self.QA, sA)
        kB = self.cB + self.linB @ sB + np.einsum('i,cij,j->c', sB, self.QB, sB)
        return kA + kB

def int_left_kernel(Mat):
    """integer left kernel basis of Mat (r x m): rows w with w @ Mat = 0, via fraction-free row reduction of [Mat | I] (Python ints, HNF-like)."""
    r, m = Mat.shape
    rows = [[int(x) for x in Mat[i]] + [1 if j == i else 0 for j in range(r)] for i in range(r)]
    piv_rows = []
    for col in range(m):
        # gcd-reduce column col among remaining rows
        while True:
            nz = [i for i in range(len(rows)) if rows[i][col] != 0]
            if len(nz) <= 1: break
            nz.sort(key=lambda i: abs(rows[i][col]))
            p = nz[0]
            for i in nz[1:]:
                q = rows[i][col] // rows[p][col]
                rows[i] = [x - q * y for x, y in zip(rows[i], rows[p])]
        nz = [i for i in range(len(rows)) if rows[i][col] != 0]
        if nz:
            piv_rows.append(rows.pop(nz[0]))
    # remaining rows have zero Mat part
    return [r_[m:] for r_ in rows]

def linear_system(P, T=None):
    """P: Poly.  Returns (E_A (e x nA), E_B (e x nB), rhs (e,), W (e x 27), const_rows) so that every full solution of key = T satisfies
       E_A sA + E_B sB = rhs; plus 'consistency' rows (kernel combos with zero linear part) with their residuals (must all be 0 for a solution to exist)."""
    D = P.S.D
    if T is None: T = np.full(D, -4, dtype=np.int64)
    Qall = np.concatenate([P.QA.reshape(D, -1), P.QB.reshape(D, -1)], axis=1)
    keep = np.flatnonzero(np.any(Qall != 0, axis=0)); Qall = Qall[:, keep]
    K = int_left_kernel(Qall)                      # integer lattice basis of combos w with w.Qall = 0
    K = np.array(K, dtype=object) if K else np.zeros((0, D), dtype=object)
    lin = np.concatenate([P.linA, P.linB], axis=1).astype(object)
    const = (P.cA + P.cB).astype(object)
    rows = []; res = []
    for w in K:
        l = w @ lin; rhs = w @ (T.astype(object) - const)
        rows.append((l, rhs))
    return K, rows

if __name__ == '__main__':
    S = build(3); LQ = make(S, 'Z2^2xZ13')
    seed = sys.argv[1] if len(sys.argv) > 1 else '3'
    w = np.load('wit_Z22xZ13_u3_s%s.npy' % seed); a, b = w[:LQ.RK], w[LQ.RK:]
    P = Poly(S, LQ, a, b)
    print('nA', P.nA, 'nB', P.nB)
    rng = np.random.default_rng(5)
    for t in range(200):
        sA = rng.choice([-1, 1], P.nA); sB = rng.choice([-1, 1], P.nB)
        ke, XA, XB = P.keys_exact(sA, sB); kp = P.keys_poly(sA, sB)
        assert np.array_equal(ke, kp), (t, ke, kp)
    print('poly decomposition matches exact keys on 200 random signs (all 27 classes)')
    K, rows = linear_system(P)
    print('kernel dim', len(K))
    lin_rank = np.linalg.matrix_rank(np.array([np.array(l, dtype=float) for l, r in rows])) if rows else 0
    print('linear rank (A+B joint)', lin_rank)

from math import gcd
def reduce_linear(P, T=None):
    """Integer-reduced linear system E s = rhs (s = (sA, sB)), e equations; returns (E (e x n) int64, rhs (e,) int64, ok, info).  ok False => no lift can exist
       (inconsistent or parity-impossible).  Row ops are unimodular over Z on the kernel-lattice rows, then rows are divided by the gcd of their lin entries."""
    K, rows = linear_system(P, T)
    n = P.nA + P.nB
    Aug = [[int(x) for x in l] + [int(r)] for l, r in rows]
    # integer echelon on the lin columns
    piv = []
    for col in range(n):
        while True:
            nz = [i for i in range(len(Aug)) if Aug[i][col] != 0 and i not in [p for p in piv]]
            if len(nz) <= 1: break
            nz.sort(key=lambda i: abs(Aug[i][col])); p = nz[0]
            for i in nz[1:]:
                q = Aug[i][col] // Aug[p][col]
                Aug[i] = [x - q * y for x, y in zip(Aug[i], Aug[p])]
        nz = [i for i in range(len(Aug)) if Aug[i][col] != 0 and i not in piv]
        if nz: piv.append(nz[0])
    rest = [Aug[i] for i in range(len(Aug)) if i not in piv]
    bad = [r for r in rest if r[n] != 0]
    E = []; rhs = []; ok = not bad; par = True
    for i in piv:
        r = Aug[i]; g = 0
        for x in r[:n]: g = gcd(g, abs(x))
        if r[n] % g != 0: par = False
        E.append([x // g for x in r[:n]]); rhs.append(r[n] // g if r[n] % g == 0 else r[n] / g)
    return (np.array(E, dtype=np.int64), np.array(rhs), ok and par, dict(inconsistent=len(bad), parity_ok=par, e=len(E), kernel=len(K)))

def natural_linear(P, T=None):
    """Equations from pairs/groups of classes with IDENTICAL quadratic part: key_c - key_c' = T_c - T_c' with a purely linear left side.
       Returns (E (e x n) int64 with rows reduced by gcd, rhs (e,) int64, ok, info).  ok False if some equation is parity-impossible (no lift exists).
       Rank-completeness against the full kernel-derived rank is reported in info."""
    D = P.S.D
    if T is None: T = np.full(D, -4, dtype=np.int64)
    Qall = np.concatenate([P.QA.reshape(D, -1), P.QB.reshape(D, -1)], axis=1)
    lin = np.concatenate([P.linA, P.linB], axis=1); const = P.cA + P.cB
    groups = {}
    for c in range(D): groups.setdefault(Qall[c].tobytes(), []).append(c)
    rows = []; rhs = []; ok = True
    for key, cs in groups.items():
        for c in cs[1:]:
            l = lin[c] - lin[cs[0]]; r = (T[c] - const[c]) - (T[cs[0]] - const[cs[0]])
            if not l.any():
                if r != 0: ok = False
                continue
            g = int(np.gcd.reduce(np.abs(l)))
            if r % g != 0: ok = False
            rows.append(l // g); rhs.append(r // g if r % g == 0 else r)
    E = np.array(rows, dtype=np.int64).reshape(-1, lin.shape[1]); rhs = np.array(rhs, dtype=np.int64)
    # rank check against full kernel-derived system
    K, krows = linear_system(P, T)
    if krows:
        full = np.array([np.array(l, dtype=float) for l, r in krows]); rk_full = np.linalg.matrix_rank(full)
        aug = np.array([np.array(list(l) + [r], dtype=float) for l, r in krows]); rk_aug = np.linalg.matrix_rank(aug)
    else: rk_full = rk_aug = 0
    rkE = np.linalg.matrix_rank(E.astype(float)) if len(E) else 0
    if rk_aug != rk_full: ok = False
    # dependent natural rows would be redundant: drop to an independent subset (first rkE independent rows)
    return E, rhs, ok, dict(n_nat=len(E), rank_nat=int(rkE), rank_full=int(rk_full), rank_aug=int(rk_aug), groups=sorted(len(v) for v in groups.values()))

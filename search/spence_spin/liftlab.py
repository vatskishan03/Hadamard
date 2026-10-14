"""Lift a Q = Z_2^2 x Z_13 level solution (a, b) to sign vectors on G = Z_2^2 x Z_11 x Z_13 (fibre Z_11, H = QR_11).  Q-values {+-1}: 'Legendre' fibre (b0, s, -s), free sign s;
{+-9, +-11}: 'spike' fibre, determined.  Provides exact lift-lattice test (full 27-class keys over sampled lifts)."""
import numpy as np, sys, time
from struct573 import build
from levels import make
from lattice import Lat

QR = {1, 3, 4, 5, 9}
def triple(val, s):
    return {1: (1, s, -s), -1: (-1, s, -s), 9: (-1, 1, 1), 11: (1, 1, 1), -9: (1, -1, -1), -11: (-1, -1, -1)}[int(val)]

class Lifter:
    def __init__(self, S, LQ):
        self.S, self.LQ = S, LQ
        G = S.G
        self.clsG = np.zeros(G.v, dtype=np.int64)
        r = G.coords[:, 2] % 11
        self.clsG = np.where(r == 0, 0, np.where(np.isin(r, list(QR)), 1, 2))
        self.repB = [int(np.flatnonzero(S.orbH == o)[0]) for o in range(S.kH)]
        self.repA = [int(np.flatnonzero(S.orbK == o)[0]) for o in range(S.kK)]
        self.qB = np.array([LQ.cos[g] for g in self.repB]); self.cB = self.clsG[self.repB]
        self.qA = np.array([LQ.orbK[LQ.cos[g]] for g in self.repA]); self.cA = self.clsG[self.repA]
    def legendre_points(self, a, b):
        return [i for i in range(len(b)) if abs(b[i]) == 1], [R for R in range(len(a)) if abs(a[R]) == 1]
    def lift(self, a, b, sA, sB):
        """sA: dict Q-orbit -> +-1 ; sB: dict Q-point -> +-1 (only for Legendre points)"""
        XB = np.zeros(self.S.kH); XA = np.zeros(self.S.kK)
        tb = {q: triple(b[q], sB.get(q, 1)) for q in range(len(b))}
        ta = {R: triple(a[R], sA.get(R, 1)) for R in range(len(a))}
        for o in range(self.S.kH): XB[o] = tb[int(self.qB[o])][int(self.cB[o])]
        for o in range(self.S.kK): XA[o] = ta[int(self.qA[o])][int(self.cA[o])]
        return XA, XB
    def check_compress(self, a, b, XA, XB):
        S, LQ = self.S, self.LQ
        cntA = np.zeros((S.kK, LQ.nQ)); np.add.at(cntA, (S.orbK, LQ.cos), 1)
        cntB = np.zeros((S.kH, LQ.nQ)); np.add.at(cntB, (S.orbH, LQ.cos), 1)
        repA = [int(np.flatnonzero(LQ.orbK == r)[0]) for r in range(LQ.RK)]; repB = [int(np.flatnonzero(LQ.orbH == r)[0]) for r in range(LQ.RH)]
        return np.array_equal(np.rint(XA @ cntA[:, repA]), a) and np.array_equal(np.rint(XB @ cntB[:, repB]), b)

def lattice_test(S, LQ, a, b, nsamp=3000, seed=1):
    lf = Lifter(S, LQ); rng = np.random.default_rng(seed)
    LB, LA = lf.legendre_points(a, b)[0], lf.legendre_points(a, b)[1]
    def rand_s(): return ({R: int(rng.choice([-1, 1])) for R in LA}, {q: int(rng.choice([-1, 1])) for q in LB})
    sA0, sB0 = rand_s(); XA0, XB0 = lf.lift(a, b, sA0, sB0)
    assert lf.check_compress(a, b, XA0, XB0)
    KA0 = np.rint(S.valsA(XA0[None])).astype(np.int64)[0]; KB0 = np.rint(S.valsB(XB0[None])).astype(np.int64)[0]
    La, Lb = Lat(S.D), Lat(S.D)
    KAs, KBs = [], []
    for i in range(nsamp):
        sA, sB = rand_s(); XA, XB = lf.lift(a, b, sA, sB)
        KAs.append(XA); KBs.append(XB)
    KA = np.rint(S.valsA(np.array(KAs))).astype(np.int64); KB = np.rint(S.valsB(np.array(KBs))).astype(np.int64)
    # A and B lifts are independent: lattice of A-differences and B-differences
    for i in range(1, nsamp): La.add(KA[i] - KA[0]); Lb.add(KB[i] - KB[0])
    Lab = Lat(S.D)
    for L in (La, Lb):
        for c, r in L.rows.items(): Lab.add(r)
    # key of base lift
    T0 = KA[0] + KB[0]
    tgt = np.full(S.D, -4) - T0
    return dict(nLegendreA=len(LA), nLegendreB=len(LB), rank=len(Lab.rows), diag=Lab.index_data(), target_in_lattice=Lab.contains(tgt),
                residual_class_range=(int(tgt.min()), int(tgt.max())), n_nonzero_resid=int((tgt != 0).sum()))

if __name__ == '__main__':
    S = build(3); LQ = make(S, 'Z2^2xZ13')
    w = np.load('wit_Z22xZ13_u3.npy'); a, b = w[:LQ.RK], w[LQ.RK:]
    t = time.time()
    print(lattice_test(S, LQ, a, b), '%.1fs' % (time.time() - t))

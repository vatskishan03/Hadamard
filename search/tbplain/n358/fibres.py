#!/usr/bin/env python3
"""Fibre utilities for G = Z2^4 x Z11: letters, cell maps, the complete A/B fibres of an R-tuple, a rejection sampler for C/D.

A letter is the vector of the six class signs of one Z11 coordinate z (64 letters, +-1 on the six E-classes of sizes
1, 3, 3, 3, 3, 3); a C- or D-sequence is a word of 11 letters, an A- or B-sequence a word of 6 letters (z and -z carry the
same letter).  The weight of a letter is w = sum_c |P_c| x_c; the vector w(z) of the column weights is the next
compression of the tower (the "Z11 level").  Cells: C has 66 cells (class c, z in Z11), A has 36 cells (class c, |z|).
`Fib.enum_A` lists every A-type sequence with a given R-vector (the complete fibre); `Fib.sample_CD` draws C/D members
with given row sums r and column weights w by rejection.  Running the file checks both against the sequence conventions.

    python3 fibres.py
"""
import itertools, math
import numpy as np
from structure import Struct
from rlevel import RLevel

WEIGHTS = np.array([1, 3, 3, 3, 3, 3])
LET = np.array(list(itertools.product([1, -1], repeat=6)), dtype=np.int8)        # 64 letters x 6 classes
LETW = LET.astype(np.int64) @ WEIGHTS                                            # weight of each letter
WVALS = sorted(set(LETW.tolist()))

class Fib:
    def __init__(self, S, RL):
        self.S, self.RL = S, RL
        p = int(S.mods[-1]); self.p = p
        self.cellC = np.full((6, p), -1, dtype=np.int64); self.cellA = np.full((6, p // 2 + 1), -1, dtype=np.int64)
        G = S.coords
        for g in range(S.m):
            c = int(RL.cls["C"][RL.piE[g]]); z = int(G[g, -1])
            self.cellC[c, z] = S.orbC[g]
            ca = int(RL.cls["A"][RL.piE[g]]); za = min(z, p - z)
            self.cellA[ca, za] = S.orbA[g]
        assert (self.cellC >= 0).all() and (self.cellA >= 0).all()
        assert len(set(self.cellC.ravel().tolist())) == 66 and len(set(self.cellA.ravel().tolist())) == S.nA
        self.byw = {w: np.flatnonzero(LETW == w) for w in WVALS}
        self.nw = {w: len(v) for w, v in self.byw.items()}
        # class sizes (E-classes) for the weights: check WEIGHTS == |P_c|
        self.csize = np.bincount(RL.cls["C"])
        assert sorted(self.csize.tolist()) == [1, 3, 3, 3, 3, 3]

    # ----- cell signs from x[c][z]
    def cells_C(self, X):          # X (n, 6, p) -> (n, 66)
        out = np.zeros((len(X), 66), dtype=np.int8)
        for c in range(6):
            for z in range(self.p): out[:, self.cellC[c, z]] = X[:, c, z]
        return out
    def cells_A(self, X):          # X (n, 6, 6)
        out = np.zeros((len(X), self.S.nA), dtype=np.int8)
        for c in range(6):
            for z in range(self.p // 2 + 1): out[:, self.cellA[c, z]] = X[:, c, z]
        return out
    # ----- A,B: complete enumeration of the R-fibre with row sums r
    def enum_A(self, r):
        rows = []
        S6 = np.array(list(itertools.product([1, -1], repeat=6)), dtype=np.int8)        # (x(0), x(1..5))
        val = S6[:, 0].astype(np.int64) + 2 * S6[:, 1:].sum(1)
        for c in range(6):
            rows.append(S6[val == r[c]])
        sizes = [len(x) for x in rows]
        n = int(np.prod(sizes))
        idx = np.indices(sizes).reshape(6, -1)
        X = np.stack([rows[c][idx[c]] for c in range(6)], 1)                              # (n, 6, 6)
        return X
    # ----- (r,w)-conditioned rejection sampler for C/D
    def sample_CD(self, r, w, n, rng, max_tries=200):
        """n members with row sums r (6) and column weights w (p); returns X (n,6,p) and acceptance rate"""
        p = self.p
        out = []; got = 0; tries = 0; acc_n = 0; tot_n = 0
        chunk = max(2000, min(500000, n * 20))
        while got < n and tries < max_tries:
            idx = np.zeros((chunk, p), dtype=np.int64)
            for z in range(p):
                bw = self.byw[int(w[z])]
                idx[:, z] = bw[rng.integers(0, len(bw), size=chunk)]
            X = LET[idx].transpose(0, 2, 1)                                               # (chunk, 6, p)
            ok = (X.sum(2) == np.asarray(r)[None, :]).all(1)
            out.append(X[ok]); got += int(ok.sum()); acc_n += int(ok.sum()); tot_n += chunk; tries += 1
        X = np.concatenate(out)[:n] if out else np.zeros((0, 6, p), dtype=np.int8)
        return X, acc_n / max(tot_n, 1)
    def wlog_n(self, w):
        return float(sum(math.log10(self.nw[int(x)]) for x in w))
    def wvec_A(self, X):            # X (n,6,6): w for z=0..5 (symmetric)
        return np.einsum("ncz,c->nz", X.astype(np.int64), WEIGHTS)
    def wvec_C(self, X):
        return np.einsum("ncz,c->nz", X.astype(np.int64), WEIGHTS)

if __name__ == "__main__":
    S = Struct("i"); RL = RLevel(S); F = Fib(S, RL)
    print("n letters per weight:", F.nw)
    rng = np.random.default_rng(1)
    # check: enum_A gives the complete fibre of an R-vector (member count = product of the per-class counts)
    r = [-1, -1, 1, -1, 3, 3]
    X = F.enum_A(r)
    assert ((X[:, :, 0] + 2 * X[:, :, 1:].sum(2)) == np.asarray(r)[None, :]).all() and len(np.unique(X.reshape(len(X), -1), axis=0)) == len(X)
    print("enum_A ok: %d members for r = %s" % (len(X), r))
    # check: sampler output is in the R-fibre and has the required w
    Xc = rng.choice([-1, 1], size=(1, 6, 11)).astype(np.int8)
    r2 = Xc.sum(2)[0]; w2 = F.wvec_C(Xc)[0]
    Y, acc = F.sample_CD(r2, w2, 2000, rng)
    assert (Y.sum(2) == r2).all() and (F.wvec_C(Y) == w2).all()
    print("sample_CD ok: acceptance %.4g, log10(prod n(w_z)) %.2f -> fibre estimate 10^%.2f" % (acc, F.wlog_n(w2), F.wlog_n(w2) + math.log10(acc)))
    # check cells -> sequences: class sums
    cs = F.cells_C(Y[:5]); full = S.expand(cs, "C")
    # compare R-vector computed by RLevel conventions: row sums over z of class values must equal compression of the sequence over Z11
    comp = full.reshape(5, 16, 11).sum(2)                                           # (5,16) over E index
    cl = RL.cls["C"]
    assert all((comp[:, cl == c] == r2[c]).all() for c in range(6)), "R-vector mismatch"
    print("cells/sequence consistency ok")

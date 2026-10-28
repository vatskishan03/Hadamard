#!/usr/bin/env python3
"""The structure behind QH(358): G = Z2^4 x Z11, H = <block Singer> of order 3, Q = Z2^2, T = the three nontrivial characters.

Written from complex-hadamard/TBORDERED.md.  G = Z_m1 x ... x Z_mr (mixed radix, last factor fastest), H = <M> with M an
integer matrix acting on coordinate column vectors, phi: G -> Z2^2 given by an integer matrix Lphi (rows mod 2).  Four
+-1 sequences A, B, C, D on G with A, B symmetric and constant on <H, -1>-orbits and C, D constant on H-orbits satisfy
sum_k PAF_k(g) = -Gamma(phi(g)) for g != 0 and 4m at g = 0, Gamma = 12 on ker(phi) and -4 elsewhere.  `Struct` holds the
orbits, the PAF classes (one per <H, -1>-orbit), the target and the exact PAF (FFT, rounded, or by the definition).
Structures: "i" = Z2^4 x Z11 (N = 179, order 358), "i3" = Z2^4 x Z3 (N = 51, order 102; the small analogue of
small_analogue.py).  Running the file checks the class tables and the FFT PAF against the direct PAF on random sequences.

    python3 structure.py
"""
import itertools
import numpy as np

STRUCTS = {
    "i":  dict(mods=[2, 2, 2, 2, 11], M=[[0,1,0,0,0],[1,1,0,0,0],[0,0,0,1,0],[0,0,1,1,0],[0,0,0,0,1]], Lphi=[[0,0,1,0,0],[0,0,0,1,0]]),
    "i3": dict(mods=[2, 2, 2, 2, 3], M=[[0,1,0,0,0],[1,1,0,0,0],[0,0,0,1,0],[0,0,1,1,0],[0,0,0,0,1]], Lphi=[[0,0,1,0,0],[0,0,0,1,0]]),
}


class Struct:
    def __init__(self, name):
        d = STRUCTS[name]
        self.name = name
        self.mods = np.array(d["mods"]); self.r = len(self.mods)
        self.M = np.array(d["M"], dtype=np.int64); self.Lphi = np.array(d["Lphi"], dtype=np.int64)
        self.m = int(np.prod(self.mods)); self.t = 3; self.N = self.m + self.t
        self.shape = tuple(int(x) for x in self.mods)
        self.coords = np.array(list(itertools.product(*[range(int(x)) for x in self.mods])), dtype=np.int64)
        self.idx = lambda c: np.ravel_multi_index(tuple(((c[..., i] % self.mods[i])) for i in range(self.r)), self.shape)
        assert (self.idx(self.coords) == np.arange(self.m)).all()
        self.hperm = self.idx(self.coords @ self.M.T)
        assert len(set(self.hperm.tolist())) == self.m
        # order of H
        p = np.arange(self.m); o = 0
        while True:
            p = self.hperm[p]; o += 1
            if (p == np.arange(self.m)).all(): break
        self.order = o
        self.neg = self.idx(-self.coords)
        self.orbC = self._orbits([self.hperm])
        self.orbA = self._orbits([self.hperm, self.neg])
        self.nC, self.nA = int(self.orbC.max()) + 1, int(self.orbA.max()) + 1
        self.sizeC = np.bincount(self.orbC); self.sizeA = np.bincount(self.orbA)
        self.phi = (self.coords @ self.Lphi.T) % 2                      # (m, 2)
        self.q = self.phi[:, 0] * 2 + self.phi[:, 1]                    # Q index 0..3
        gam = np.array([12, -4, -4, -4])
        self.target = -gam[self.q]; self.target[0] = 4 * self.m
        self.repA = np.array([np.flatnonzero(self.orbA == c)[0] for c in range(self.nA)])
        self.tvecA = self.target[self.repA]
        for c in range(self.nA):
            assert (self.target[self.orbA == c] == self.tvecA[c]).all(), "target not constant on a class"
        self.zero_class = int(self.orbA[0])
        self.cols = [c for c in range(self.nA) if c != self.zero_class]    # the 35 (or 17) nonzero classes
        self.T = self.tvecA[self.cols].astype(np.int64)

    def _orbits(self, perms):
        m = self.m; par = list(range(m))
        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]; x = par[x]
            return x
        for pm in perms:
            for i in range(m):
                a, b = find(i), find(int(pm[i]))
                if a != b: par[a] = b
        lab = {}; out = np.zeros(m, dtype=np.int64)
        for i in range(m):
            r = find(i)
            out[i] = lab.setdefault(r, len(lab))
        return out

    # ---- exact PAF
    def paf_direct(self, X):
        """PAF over G by the definition (slow, exact, integer)"""
        X = np.asarray(X, dtype=np.int64)
        add = self.idx(self.coords[:, None, :] + self.coords[None, :, :])      # add[g,h] = g + h
        return np.array([int((X * X[add[g]]).sum()) for g in range(self.m)])

    def paf_classes(self, S):
        """S: (n, m) int array of +-1 -> (n, nA) integer PAF values at the class representatives (FFT, rounded)"""
        S = np.asarray(S)
        n = S.shape[0]; out = np.zeros((n, self.nA), dtype=np.int32)
        for s in range(0, n, 20000):
            blk = S[s:s + 20000].reshape((-1,) + self.shape).astype(np.float64)
            F = np.fft.fftn(blk, axes=tuple(range(1, 1 + self.r)))
            P = np.fft.ifftn(np.abs(F) ** 2, axes=tuple(range(1, 1 + self.r))).real.reshape(blk.shape[0], -1)
            out[s:s + 20000] = np.rint(P[:, self.repA]).astype(np.int32)
        return out

    def expand(self, s, kind):
        """cell signs (n, ncells) -> sequences over G (n, m)"""
        return s[:, self.orbA if kind == "A" else self.orbC]

    def describe(self):
        return dict(name=self.name, mods=self.mods.tolist(), m=self.m, N=self.N, order_H=self.order, nC=self.nC, nA=self.nA,
                    sizesC=sorted(set(self.sizeC.tolist())), sizesA=sorted(set(self.sizeA.tolist())), ncols=len(self.cols))


if __name__ == "__main__":
    for nm in ("i", "i3"):
        S = Struct(nm)
        print(S.describe())
        print(" class counts: nA", S.nA, "nC", S.nC, " target classes (nonzero) values", sorted(set(S.T.tolist())), " zero class target", S.tvecA[S.zero_class])
        # exact PAF check on random sequences
        rng = np.random.default_rng(1)
        sA = rng.choice([-1, 1], size=(3, S.nA)); sC = rng.choice([-1, 1], size=(3, S.nC))
        for kind, s in (("A", sA), ("C", sC)):
            X = S.expand(s, kind)
            F = S.paf_classes(X)
            for i in range(3):
                d = S.paf_direct(X[i])
                assert (d[S.repA] == F[i]).all(), "FFT PAF != direct PAF"
                # PAF constant on classes (H-invariant + symmetric for A; for C only H-orbit constancy + PAF(g)=PAF(-g))
                for c in range(S.nA):
                    assert len(set(d[S.orbA == c].tolist())) == 1
                if kind == "A":
                    assert (X[i] == X[i][S.neg]).all()
        print(" PAF FFT == direct convolution, class-constant: OK")

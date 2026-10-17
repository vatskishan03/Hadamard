"""Spence-spin core: G abelian of order v = n-1, H <= Aut(G), mu in Aut(G) with mu normalising H and mu^3 acting inside H-orbits.
Quadruple (A, B, B o mu, B o mu^2):  A is K=<H,mu>-invariant with row sum 2, B is H-invariant with row sum 0,
    PAF_A(g) + PAF_B(g) + PAF_B(mu g) + PAF_B(mu^2 g) = -4   for g != 0         (Spence, t = 1)
The second pattern (A, A o mu, A o mu^2, D) is the same problem with the blocks renamed (the Spence condition and the
row sums (2,0,0,0) are symmetric in the four blocks), so only the first is implemented.
"""
import math
import numpy as np
from math import gcd, log10, pi


class Grp:
    def __init__(self, mods):
        self.mods = np.array(mods, dtype=np.int64)
        self.r = len(mods)
        self.v = int(np.prod(self.mods))
        self.coords = np.array(np.unravel_index(np.arange(self.v), tuple(int(m) for m in mods))).T.astype(np.int64).reshape(self.v, self.r)
        self.w = np.array([int(np.prod(self.mods[i + 1:])) for i in range(self.r)], dtype=np.int64)

    def index(self, C):
        return (np.asarray(C, dtype=np.int64) % self.mods) @ self.w

    def neg(self):
        return self.index(-self.coords)

    def aut(self, M):
        """x -> M x ; M integer r x r (column j = image of generator e_j).  Checked: well defined, bijective."""
        M = np.asarray(M, dtype=np.int64).reshape(self.r, self.r)
        if np.any((self.mods[None, :] * M) % self.mods[:, None] != 0):
            raise ValueError('matrix not well defined on G')
        p = self.index(self.coords @ M.T)
        if len(np.unique(p)) != self.v:
            raise ValueError('not bijective')
        return p

    def shift(self, s):
        return self.index(self.coords + self.coords[s])


def closure_partition(v, perms):
    lab = -np.ones(v, dtype=np.int64)
    k = 0
    for s in range(v):
        if lab[s] >= 0:
            continue
        lab[s] = k
        st = [s]
        while st:
            x = st.pop()
            for p in perms:
                y = int(p[x])
                if lab[y] < 0:
                    lab[y] = k
                    st.append(y)
        k += 1
    return lab, k


class Spin:
    """The structure (G, H, mu)."""

    def __init__(self, G, Hperms, mu, name=''):
        self.G, self.name = G, name
        v = G.v
        self.Hperms, self.mu = list(Hperms), mu
        self.mu2 = mu[mu]
        self.mu3 = mu[self.mu2]
        self.orbH, self.kH = closure_partition(v, self.Hperms)
        self.orbH3, k3 = closure_partition(v, self.Hperms + [self.mu3])
        if k3 != self.kH:
            raise ValueError('mu^3 does not preserve the H-orbits')
        # mu maps H-orbits into H-orbits
        img = np.array([self.orbH[self.mu[np.flatnonzero(self.orbH == o)[0]]] for o in range(self.kH)])
        if not np.array_equal(img[self.orbH], self.orbH[self.mu]):
            raise ValueError('mu does not normalise H')
        self.orbK, self.kK = closure_partition(v, self.Hperms + [mu])
        self.trivial_spin = (self.kK == self.kH)          # mu acts inside H-orbits: not a real spin structure
        self.szH = np.bincount(self.orbH, minlength=self.kH)
        self.szK = np.bincount(self.orbK, minlength=self.kK)
        self.neg = G.neg()
        cl, _ = closure_partition(v, self.Hperms + [mu, self.neg])
        z = cl[0]
        labs = sorted(set(cl.tolist()) - {z})
        remap = {l: i for i, l in enumerate(labs)}
        self.cls = np.array([remap.get(int(c), -1) for c in cl])
        self.D = len(labs)
        self.reps = [int(np.flatnonzero(self.cls == c)[0]) for c in range(self.D)]
        self.csize = np.bincount(self.cls[self.cls >= 0], minlength=self.D)
        self._MA = self._MB = None

    # ---- class matrices -------------------------------------------------------
    def _paf_matrix(self, orb, k, g):
        M = np.zeros((k, k), dtype=np.int64)
        np.add.at(M, (orb, orb[self.G.shift(g)]), 1)
        return M

    @property
    def MA(self):
        if self._MA is None:
            self._MA = np.array([self._paf_matrix(self.orbK, self.kK, g) for g in self.reps])
        return self._MA

    @property
    def MB(self):
        if self._MB is None:
            out = []
            for g in self.reps:
                M = self._paf_matrix(self.orbH, self.kH, g) + self._paf_matrix(self.orbH, self.kH, int(self.mu[g])) \
                    + self._paf_matrix(self.orbH, self.kH, int(self.mu2[g]))
                out.append(M)
            self._MB = np.array(out)
        return self._MB

    def valsA(self, X):
        return np.stack([np.einsum('ni,ij,nj->n', X, self.MA[c], X) for c in range(self.D)], axis=1)

    def valsB(self, X):
        return np.stack([np.einsum('ni,ij,nj->n', X, self.MB[c], X) for c in range(self.D)], axis=1)

    # ---- row sums ---------------------------------------------------------------
    @staticmethod
    def _dp(sizes, v):
        O = len(sizes)
        W = np.zeros((O + 1, 2 * v + 1))
        W[O, v] = 1.0
        for j in range(O - 1, -1, -1):
            s = int(sizes[j])
            W[j, s:] += W[j + 1, :2 * v + 1 - s]
            W[j, :2 * v + 1 - s] += W[j + 1, s:]
        return W

    def pool(self, which, r):
        sizes = self.szK if which == 'A' else self.szH
        W = self._dp(sizes, self.G.v)
        return W[0, r + self.G.v]

    def sample(self, which, r, count, rng):
        sizes = self.szK if which == 'A' else self.szH
        v = self.G.v
        W = self._dp(sizes, v)
        if W[0, r + v] == 0:
            return None
        O = len(sizes)
        Y = np.empty((count, O))
        rem = np.full(count, r)
        for j in range(O):
            s = int(sizes[j])
            tot = W[j, rem + v]
            plus = np.where(rem - s >= -v, W[j + 1, np.clip(rem - s + v, 0, 2 * v)], 0.0)
            take = rng.random(count) < plus / tot
            Y[:, j] = np.where(take, 1.0, -1.0)
            rem = rem - Y[:, j].astype(np.int64) * s
        assert np.all(rem == 0)
        return Y

    # ---- Gaussian estimates --------------------------------------------------------
    def estimate(self, samples=1500, seed=1):
        """log10 of the expected number of solutions (A,B), Gaussian local-limit model: per-class diagonal and full covariance.
        Both ignore correlations between the lattice conditions; unreliable for composite G."""
        NA, NB = self.pool('A', 2), self.pool('B', 0)
        if NA == 0 or NB == 0:
            return dict(lE=-math.inf, lEfull=-math.inf, NA=NA, NB=NB)
        rng = np.random.default_rng(seed)
        XA = self.sample('A', 2, samples, rng)
        XB = self.sample('B', 0, samples, rng)
        VA = np.rint(self.valsA(XA)).astype(np.int64)
        VB = np.rint(self.valsB(XB)).astype(np.int64)
        # lattice: all values VA+VB[c] are congruent to (VA0 + VB0) mod g
        lat = []
        for c in range(self.D):
            g = 0
            for V in (VA, VB):
                u = np.unique(V[:, c])
                if len(u) > 1:
                    g = gcd(g, int(np.gcd.reduce(np.diff(u))))
            r0 = int(VA[0, c] + VB[0, c])
            lat.append((g, r0))
        mu_ = VA.mean(0) + VB.mean(0)
        cov = np.cov(VA.T) + np.cov(VB.T) if self.D > 1 else np.array([[VA.var() + VB.var()]])
        cov = np.atleast_2d(cov)
        var = np.diag(cov).copy()
        tgt = -4.0
        # lattice check
        for c, (g, r0) in enumerate(lat):
            if g and (int(tgt) - r0) % g:
                return dict(lE=-math.inf, lEfull=-math.inf, NA=NA, NB=NB, lattice_fail=c)
        # drop the class carrying the weighted largest spread (row sums determine it)
        w = self.csize * np.sqrt(np.maximum(var, 1e-12))
        drop = int(np.argmax(w))
        keep = [c for c in range(self.D) if c != drop and var[c] > 1e-9]
        lp = 0.0
        for c in keep:
            g = lat[c][0] or 1
            lp += -0.5 * log10(2 * pi * var[c]) - (tgt - mu_[c]) ** 2 / (2 * var[c]) / math.log(10) + log10(g)
        lE = log10(NA) + log10(NB) + lp
        # full covariance
        lEf = None
        if keep:
            Cm = cov[np.ix_(keep, keep)]
            try:
                sign, logdet = np.linalg.slogdet(Cm)
                d = (tgt - mu_[keep])
                maha = float(d @ np.linalg.solve(Cm, d)) if sign > 0 else float('nan')
                if sign > 0:
                    lg = sum(log10(lat[c][0] or 1) for c in keep)   # diagonal lattice volume (lower bound on the true index)
                    lEf = log10(NA) + log10(NB) - 0.5 * len(keep) * log10(2 * pi) - 0.5 * logdet / math.log(10) - 0.5 * maha / math.log(10) + lg
            except Exception:
                lEf = None
        return dict(lE=lE, lEfull=lEf, NA=NA, NB=NB, lNA=log10(NA), lNB=log10(NB), keep=len(keep))

    # ---- exhaustive solver (small structures) ---------------------------------------
    def pool_vectors(self, which, r, limit=1 << 24):
        sizes = self.szK if which == 'A' else self.szH
        k = len(sizes)
        if k > 26:
            raise MemoryError('pool too large')
        # enumerate by meet in the middle on the row sum
        n = 1 << k
        if n > limit * 8:
            raise MemoryError('pool too large')
        idx = np.arange(n, dtype=np.int64)
        bits = ((idx[:, None] >> np.arange(k)[None, :]) & 1)
        X = (1 - 2 * bits).astype(np.int8)
        rs = X.astype(np.int64) @ sizes
        return X[rs == r]

    def solve_exhaustive(self, maxsol=50):
        XA = self.pool_vectors('A', 2).astype(np.float64)
        XB = self.pool_vectors('B', 0).astype(np.float64)
        VA = np.rint(np.stack([np.einsum('ni,ij,nj->n', XA, self.MA[c].astype(np.float64), XA) for c in range(self.D)], 1)).astype(np.int64)
        VB = np.rint(np.stack([np.einsum('ni,ij,nj->n', XB, self.MB[c].astype(np.float64), XB) for c in range(self.D)], 1)).astype(np.int64)
        d = {}
        for i, row in enumerate(VB):
            d.setdefault(row.tobytes(), []).append(i)
        sols = []
        tgt = -4 - VA
        for i, row in enumerate(tgt):
            hit = d.get(row.tobytes())
            if hit:
                for j in hit:
                    sols.append((i, j))
                    if len(sols) >= maxsol:
                        return sols, XA, XB
        return sols, XA, XB

    def sequences(self, xa, xb):
        """expand orbit sign vectors to full sequences (A, B, B o mu, B o mu^2) as +-1 arrays over G"""
        A = xa.astype(np.int64)[self.orbK]
        B = xb.astype(np.int64)[self.orbH]
        return [A, B, B[self.mu], B[self.mu2]]


def paf_sum_check(G, seqs):
    v = G.v
    tot = np.zeros(v, dtype=np.int64)
    for s in range(v):
        sh = G.shift(s)
        tot[s] = sum(int(x @ x[sh]) for x in seqs)
    return tot

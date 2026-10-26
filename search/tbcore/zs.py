#!/usr/bin/env python3
"""Zero-diagonal Szollosi array (plain, bordered, character-bordered) over G = prod Z_{m_i} with multiplier group H.

A, B symmetric (constant on <H,-1>-orbits) with A(0) = B(0) = 0; C, D +-1 constant on H-orbits.
sum_k PAF_k(g) = Tg(g) for g != 0, with Tg(g) = -Gamma(phi(g)), Gamma(a) = 4 sum_{chi in T} chi(a) (T closed under
conjugation); Tg == -4 for the bordered array (Q trivial, T = {1}); sum_k PAF_k(0) = 4|G| - 2.
Also general-GS mode (all four blocks constant on H-orbits, no zero at 0): kind 'GS'.
"""
import itertools, math, sys
from math import gcd, log10, pi
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


def hnf_basis(vecs, dim):
    B, piv = [], []
    for v in vecs:
        v = [int(x) for x in v]
        while True:
            j = next((c for c in range(dim) if v[c] != 0), None)
            if j is None:
                break
            if j in piv:
                k = piv.index(j); b = B[k]
                a0, a1 = b[j], v[j]
                if a1 % a0 == 0:
                    q = a1 // a0
                    v = [x - q * y for x, y in zip(v, b)]
                    continue
                g, s, t = egcd(a0, a1)
                nb = [s * x + t * y for x, y in zip(b, v)]
                nv = [(a0 // g) * y - (a1 // g) * x for x, y in zip(b, v)]
                B[k] = nb; v = nv
                continue
            if v[j] < 0:
                v = [-x for x in v]
            B.append(v); piv.append(j)
            break
    order = sorted(range(len(B)), key=lambda k: piv[k])
    return [B[k] for k in order], [piv[k] for k in order]


def egcd(a, b):
    x0, x1, y0, y1 = 1, 0, 0, 1
    while b:
        q, a, b = a // b, b, a % b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return a, x0, y0


class ZS:
    def __init__(self, mods, gens, phi=None, T=None, qmods=None, gs=False, name="", plain=False, gstarget=None, zdiag=True):
        """mods: moduli of the coordinates of G; gens: integer matrices (acting on column coordinate vectors, row i
        reduced mod mods[i]); phi: function coords(array m x d) -> labels in Z_{qmods} (array m x dq) or None;
        T: list of characters of Q given as integer tuples k (chi_k(a) = exp(2 pi i sum k_j a_j / qmods_j)).
        gs=True: four general blocks (Goethals-Seidel type), A(0) not forced to 0."""
        self.name = name
        self.mods = tuple(int(x) for x in mods)
        d = len(self.mods)
        self.d = d
        self.m = int(np.prod(self.mods))
        self.gs = gs
        self.zdiag = zdiag
        self.coords = np.array(list(itertools.product(*[range(x) for x in self.mods])), dtype=np.int64)
        mo = np.array(self.mods, dtype=np.int64)
        self.mo = mo

        def idx(c):
            return np.ravel_multi_index(tuple((c % mo).T), self.mods)

        self.idx = idx
        perms = []
        for M in gens:
            M = np.asarray(M, dtype=np.int64)
            p = idx(self.coords @ M.T)
            if len(np.unique(p)) != self.m:
                raise ValueError("generator is not an automorphism")
            perms.append(p)
        self.perms = perms
        self.neg = idx(-self.coords)
        self.orbA = self._orbits(perms)
        self.orbL = self._orbits(perms + [self.neg])
        self.nA = int(self.orbA.max()) + 1
        self.nL = int(self.orbL.max()) + 1
        # classes = L-orbits
        self.C = self.nL
        self.rep = np.array([np.flatnonzero(self.orbL == c)[0] for c in range(self.C)])
        self.cls_size = np.bincount(self.orbL)
        self.zero_class = int(self.orbL[0])
        # target function
        Tg = np.full(self.m, 4, dtype=np.int64)
        self.Gamma_info = None
        if phi is not None and T is not None:
            lab = phi(self.coords)                       # m x dq
            qm = np.array(qmods)
            gam = np.zeros(self.m)
            for k in T:
                gam += 4 * np.cos(2 * pi * ((lab * np.array(k)) / qm).sum(1))
            assert np.allclose(gam, np.rint(gam), atol=1e-6)
            Tg = -np.rint(gam).astype(np.int64)
            self.qlab = lab
            self.t = len(T)
        else:
            Tg = np.full(self.m, 0 if plain else -4, dtype=np.int64)
            self.t = 1
        Tg = Tg.copy()
        Tg[0] = (4 * self.m - 2) if (not gs and zdiag) else 4 * self.m
        if gstarget is not None:
            Tg = np.array(gstarget, dtype=np.int64)
        self.Tg = Tg
        self.tvec = np.zeros(self.C, dtype=np.int64)
        for c in range(self.C):
            v = np.unique(Tg[self.orbL == c])
            if len(v) != 1:
                raise ValueError("H does not preserve the target")
            self.tvec[c] = v[0]
        # paired / unpaired H-orbits
        self.selfpaired = 0
        self.unpaired = 0
        for o in range(self.nA):
            e = np.flatnonzero(self.orbA == o)[0]
            if self.orbA[self.neg[e]] == o:
                self.selfpaired += 1
            else:
                self.unpaired += 1

    def _orbits(self, perms):
        n = self.m
        rows = np.concatenate([np.arange(n)] * len(perms)) if perms else np.zeros(0, int)
        cols = np.concatenate(perms) if perms else np.zeros(0, int)
        g = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
        nc, lab = connected_components(g, directed=False)
        # relabel by first occurrence
        order = {}
        out = np.zeros(n, dtype=np.int64)
        for i in range(n):
            l = lab[i]
            if l not in order:
                order[l] = len(order)
            out[i] = order[l]
        return out

    # ------------------------------------------------------------------ block kinds
    def kind_orbits(self, kind):
        """returns (orbit label per element, orbit sizes, list of free orbit ids)"""
        orb = self.orbL if kind == "L" else self.orbA
        sz = np.bincount(orb)
        free = list(range(len(sz)))
        if kind == "L" and not self.gs and self.zdiag:
            free.remove(int(orb[0]))               # A(0) = 0
        return orb, sz, free

    def dp(self, kind):
        orb, sz, free = self.kind_orbits(kind)
        w = [int(sz[o]) for o in free]
        tot = sum(w)
        # f[i][x + tot] = number of sign assignments of orbits i.. with sum x
        f = [None] * (len(w) + 1)
        cur = np.zeros(2 * tot + 1)
        cur[tot] = 1.0
        f[len(w)] = cur
        for i in range(len(w) - 1, -1, -1):
            nxt = np.zeros(2 * tot + 1)
            nxt[w[i]:] += f[i + 1][:2 * tot + 1 - w[i]]
            nxt[:2 * tot + 1 - w[i]] += f[i + 1][w[i]:]
            f[i] = nxt
        return w, tot, f

    def counts(self, kind):
        w, tot, f = self.dp(kind)
        return {x - tot: f[0][x] for x in range(2 * tot + 1) if f[0][x] > 0}

    def sample(self, kind, r, K, rng):
        orb, sz, free = self.kind_orbits(kind)
        w, tot, f = self.dp(kind)
        if f[0][r + tot] <= 0:
            return None
        rs = np.full(K, r, dtype=np.int64)
        signs = np.zeros((K, len(w)), dtype=np.int8)
        for i in range(len(w)):
            fi = f[i][rs + tot]
            fp = f[i + 1][rs - w[i] + tot] if True else None
            # guard indexes
            ok = (rs - w[i] + tot >= 0) & (rs - w[i] + tot <= 2 * tot)
            fp = np.where(ok, f[i + 1][np.clip(rs - w[i] + tot, 0, 2 * tot)], 0.0)
            p = fp / fi
            u = rng.random(K) < p
            signs[:, i] = np.where(u, 1, -1)
            rs = np.where(u, rs - w[i], rs + w[i])
        return signs

    def seq(self, kind, signs):
        """full sequences (K x m) from orbit signs"""
        orb, sz, free = self.kind_orbits(kind)
        K = signs.shape[0]
        lut = np.zeros((K, len(sz)), dtype=np.int8)
        lut[:, free] = signs
        return lut[:, orb]

    def paf_classes(self, X):
        """class values of the PAF for sequences X (K x m)"""
        K = X.shape[0]
        out = np.zeros((K, self.C), dtype=np.int64)
        for s in range(0, K, 256):
            blk = X[s:s + 256].reshape((-1,) + self.mods).astype(float)
            F = np.fft.fftn(blk, axes=tuple(range(1, self.d + 1)))
            P = np.fft.ifftn(np.abs(F) ** 2, axes=tuple(range(1, self.d + 1))).real
            P = P.reshape(P.shape[0], -1)
            out[s:s + 256] = np.rint(P[:, self.rep]).astype(np.int64)
        return out

    def patterns(self, sumsq=None):
        kinds = ["L", "L", "A", "A"] if not self.gs else ["A"] * 4
        if sumsq is None:
            sumsq = 4 * self.m - 2 if not self.gs else 4 * self.m
        cnt = {k: self.counts(k) for k in ("L", "A")}
        lim = math.isqrt(sumsq)
        out = []
        Ls = sorted(r for r in cnt["L"] if 0 <= r <= lim)
        As = sorted(r for r in cnt["A"] if 0 <= r <= lim)
        if self.gs:
            for a, b, c in itertools.combinations_with_replacement(As, 3):
                d2 = sumsq - a * a - b * b - c * c
                if d2 < 0:
                    continue
                d = math.isqrt(d2)
                if d * d == d2 and d in cnt["A"] and d >= c:
                    out.append((a, b, c, d))
            return out
        for a in Ls:
            for b in Ls:
                if b < a:
                    continue
                for c in As:
                    d2 = sumsq - a * a - b * b - c * c
                    if d2 < c * c:
                        continue
                    d = math.isqrt(d2)
                    if d * d == d2 and d in cnt["A"]:
                        out.append((a, b, c, d))
        return out

    def estimate(self, rows, K=3000, seed=1, nd=400, verbose=False):
        """lattice Gaussian estimate of the number of solutions with the given positive row sums. Returns
        (log10 E lattice, log10 E product-of-spacings, log10 pool sizes, free classes)"""
        rng = np.random.default_rng(seed)
        kinds = ["L", "L", "A", "A"] if not self.gs else ["A"] * 4
        lN, mus, covs, base, diffs, pools = 0.0, [], [], None, [], []
        for kind, r in zip(kinds, rows):
            n = self.counts(kind).get(r, 0)
            if n == 0:
                return None
            lN += log10(n); pools.append(log10(n))
            S = self.sample(kind, r, K, rng)
            V = self.paf_classes(self.seq(kind, S))
            mus.append(V.mean(0)); covs.append(np.cov(V.T.astype(float), bias=True))
            base = V[0].copy() if base is None else base + V[0]
            idx = rng.integers(1, len(V), nd)
            diffs.extend(list(V[idx] - V[0]))
        mu = sum(mus); Sg = sum(covs)
        var = np.diag(Sg)
        tg = self.tvec.astype(float)
        fixed = [c for c in range(self.C) if var[c] <= 1e-9]
        for c in fixed:
            if abs(tg[c] - mu[c]) > 1e-6:
                return ("fixed-mismatch", c)
        free = [c for c in range(self.C) if var[c] > 1e-9]
        c0 = free[1:]
        if not c0:
            return (lN, lN, pools, 0)
        D = [[int(dv[c]) for c in c0] for dv in diffs]
        B, piv = hnf_basis(D, len(c0))
        if len(B) < len(c0):
            return ("rank", len(B), len(c0))
        covol = 1
        for k, b in enumerate(B):
            covol *= b[piv[k]]
        t = [int(self.tvec[c] - base[c]) for c in c0]
        for k, b in enumerate(B):
            j = piv[k]
            if t[j] % b[j]:
                return ("coset", lN, pools)
            q = t[j] // b[j]
            t = [x - q * y for x, y in zip(t, b)]
        if any(t):
            return ("coset", lN, pools)
        gprod = 1
        for ci, c in enumerate(c0):
            g = 0
            for dv in D:
                g = gcd(g, int(dv[ci]))
            gprod *= max(g, 1)
        dvec = (tg - mu)[c0]
        Sk = Sg[np.ix_(c0, c0)]
        sign, logdet = np.linalg.slogdet(Sk)
        quad = float(dvec @ np.linalg.solve(Sk, dvec))
        dens = (-0.5 * len(c0) * math.log(2 * pi) - 0.5 * logdet - 0.5 * quad) / math.log(10)
        covol = abs(covol)
        return (lN + dens + log10(covol), lN + dens + log10(gprod), pools, len(c0))

    def info(self):
        return (f"{self.name} |G|={self.m} orbits A-kind {self.nA} (self-paired {self.selfpaired}, unpaired "
                f"{self.unpaired}) L-classes {self.C}")


def cyc_units(m):
    return [u for u in range(1, m) if gcd(u, m) == 1]


def block_diag(*blocks):
    n = sum(b.shape[0] for b in blocks)
    M = np.zeros((n, n), dtype=np.int64)
    i = 0
    for b in blocks:
        k = b.shape[0]
        M[i:i + k, i:i + k] = b
        i += k
    return M


def mat(x):
    return np.array(x, dtype=np.int64)

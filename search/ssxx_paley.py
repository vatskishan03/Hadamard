#!/usr/bin/env python3
"""Exact search for complex Hadamard matrices of order 2N on Paley-product structures, N = m q.

Structure. q is a prime, q >= 5, gcd(m, q) = 1, Z_N = Z_m x Z_q, and H = {h : h = 1 (mod m), h a square (mod q)}.
With chi the quadratic character of Z_q, J the all-ones function and delta the indicator of 0, every H-invariant block
of Szöllősi's array (search/ssxx.py: A, B symmetric, C, D arbitrary, sum of PAFs 0 off 0) is

    X(a, .) = u(a) delta + w(a) J + z(a) chi,      a in Z_m,

where at each a either w(a) = 0 and u(a), z(a) = +-1 (a z-position), or z(a) = 0, w(a) = +-1 and u(a) in
{0, -2 w(a)} (a w-position). A and B are symmetric: u and w are even and z(-a) = eps z(a), eps = chi(-1). From
delta * chi = eps chi, J * chi = 0 and chi * chi = q delta - J ((f * g)(s) = sum_a f(a + s) g(a)) the condition
splits on delta, J and chi into three equations over Z_m:

    (I)   sum_k PAF_{u_k}(s) + q PAF_{z_k}(s)                 = 4N [s = 0]
    (II)  sum_k (u_k * w_k + w_k * u_k)(s) + q PAF_{w_k}(s)   = sum_k PAF_{z_k}(s)
    (III) sum_k eps (u_k * z_k)(s) + (z_k * u_k)(s)            = 0.

The compressions y_k = u_k + q w_k to Z_m satisfy sum_k PAF_{y_k}(s) = 4N [s = 0], which is (I) + q (II); given it,
(I) and (II) are equivalent. At s = 0 it reads (q + 1) N_w - 4 N_- = 4m, with N_w the number of w-positions in the
four blocks and N_- the number with u = -2w (search/NOTES.md §15). For eps = -1, position 0 of A and of B is a
w-position, so N_w >= 2.

Search (exact). The skeletons (the w-positions of the four blocks and their values y = +-q, +-(q - 2)) are listed up
to the symmetries: multipliers of Z_m, translations and reversal of C and of D, negation of a block, A <-> B and
C <-> D. Stage 1 joins, for each skeleton, the +-1 values of y at the z-positions on sum PAF_y = 0 (join4x). The
stage-1 solutions are reduced to one per orbit of the same group. Stage 2 joins, for each of them, the z_k on (II)
and (III). Every hit is lifted to Z_N and checked as a complex Hadamard matrix (ssxx.complex_hadamard). Without
--all the search stops at the first verified hit; with --all every stage-1 orbit is processed and the exact number
of stage-2 solutions is reported, so a run that ends without a hit rules the structure out.

    python3 search/ssxx_paley.py N q [--all] [--log FILE] [--threads T] [--maxsol K]
    python3 search/ssxx_paley.py N q --raw        (no symmetry reduction: counts every solution, for testing)
    python3 search/ssxx_paley.py N q --sample K [--seed S]   (K random skeletons; for structures with too many)
    python3 search/ssxx_paley.py N q --bordered ...  (the bordered array of search/bsz.py on Z_m x GF(q), m = (N - 1)/q
                                                      even: sum PAF = -4 off 0, stage-1 target -4q, stage 2 R + 4)
"""
import argparse, itertools, json, math, os, struct, subprocess, sys, tempfile, time
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import workdir as W
import ssxx


def chi_table(q):
    c = np.full(q, -1, np.int64); c[0] = 0
    for y in range(1, q):
        c[y * y % q] = 1
    return c


def gf_chi(q):
    """Quadratic character of GF(q), q = p^k, on the elements in mixed-radix order of their coordinates in the basis
    x^(k-1), ..., x, 1 (polynomials modulo a monic irreducible f, found by search). Returns (p, k, f, chi)."""
    p = next(d for d in range(2, q + 1) if q % d == 0)
    k = round(math.log(q, p))
    assert p ** k == q
    def mul(a, b, f):                               # coefficient lists, highest degree first, length k
        prod = [0] * (2 * k - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                prod[i + j] = (prod[i + j] + x * y) % p
        for i in range(len(prod) - k):              # reduce by x^k = -(f_1 x^(k-1) + ... + f_k)
            c = prod[i]
            if c:
                for j in range(1, k + 1):
                    prod[i + j] = (prod[i + j] - c * f[j - 1]) % p
                prod[i] = 0
        return prod[-k:]
    elems = [list(c) for c in itertools.product(range(p), repeat=k)]          # mixed radix, last fastest
    for f in itertools.product(range(p), repeat=k):                             # x^k + f_1 x^(k-1) + ... + f_k
        if f[-1] == 0:
            continue
        sq = {tuple(mul(e, e, f)) for e in elems[1:]}
        if len(sq) == (q - 1) // 2 and all(tuple(mul(a, b, f)) != tuple([0] * k) for a in elems[1:] for b in elems[1:]):
            chi = np.array([0] + [1 if tuple(e) in sq else -1 for e in elems[1:]], np.int64)
            return p, k, list(f), chi
    raise ValueError("no irreducible polynomial found")


def admissible(m, q, bordered=False):
    """(N_w, N_-) with (q + 1) N_w - 4 N_- = 4m (bordered: 4(m - 1)) and 0 <= N_- <= N_w."""
    out = []
    rhs = 4 * (m - 1) if bordered else 4 * m
    for nw in range(0, 4 * m + 1):
        t = (q + 1) * nw - rhs
        if t >= 0 and t % 4 == 0 and t // 4 <= nw:
            out.append((nw, t // 4))
    return out


def allpm(n):
    idx = np.arange(2 ** n, dtype=np.uint32)
    out = np.empty((2 ** n, n), np.int8)
    for j in range(n):
        out[:, j] = 1 - 2 * ((idx >> np.uint32(j)) & np.uint32(1)).astype(np.int8)
    return out


def pafs(Y, h):
    Y = Y.astype(np.int32)
    return np.stack([(Y * np.roll(Y, -s, axis=1)).sum(1) for s in range(1, h + 1)], axis=1)


def corr(F, G, s):
    """(F * G)(s) = sum_a F(a + s) G(a), rows of F against rows of G (broadcast)."""
    return (np.roll(F, -s, axis=-1) * G).sum(-1)


class Paley:
    def __init__(self, N, q, bordered=False):
        v = N - 1 if bordered else N               # bordered: the group has order N - 1 (search/bsz.py)
        assert v % q == 0 and 5 <= q < 128           # y values are stored as int8
        self.N, self.v, self.q, self.m, self.bordered = N, v, q, v // q, bordered
        self.prime = all(q % p for p in range(2, int(q ** 0.5) + 1))
        m = self.m
        assert math.gcd(m, q) == 1 and m % 2 == (0 if bordered else 1) and m >= 2
        assert self.prime or not bordered
        self.h = m // 2                              # shifts 1..h cover the nonzero shifts up to sign
        self.fixed = [0] + ([m // 2] if m % 2 == 0 else [])     # fixed points of a -> -a
        self.reps = list(range(1, (m + 1) // 2))   # representatives of the pairs {a, -a}
        self.T1 = -4 * q if bordered else 0          # stage 1: sum PAF_y(s) = T1 (s != 0)
        self.T2 = 4 if bordered else 0               # stage 2: sum PAF_z(s) = R(s) + T2
        self.eps = 1 if q % 4 == 1 else -1
        self.units = [t for t in range(1, m) if math.gcd(t, m) == 1]
        if self.prime:          # prime powers: stages 1 and 2 only depend on q and eps (no lift yet)
            self.chi = chi_table(q)
            r = next(g for g in range(2, q) if all(pow(g, (q - 1) // p, q) != 1 for p in range(2, q)
                                                    if (q - 1) % p == 0 and all(p % d for d in range(2, p))))
            self.gen = next(x for x in range(v) if x % m == 1 % m and x % q == r * r % q)
        else:
            self.gen = None
            pp, kk, self.poly, self.chi = gf_chi(q)
            self.group = [m] + [pp] * kk           # Z_m x GF(q)^+, mixed radix, last factor fastest
        a = np.arange(v)
        self.ia, self.iy = a % m, a % q            # CRT coordinates of each element of Z_v

    # ---------- skeletons ----------
    def fixed_choices(self):
        """w-sets among the fixed points of a -> -a: all of them when eps = -1 (z vanishes there), any otherwise."""
        if self.eps == -1:
            return [tuple(self.fixed)]
        return [c for r in range(len(self.fixed) + 1) for c in itertools.combinations(self.fixed, r)]

    def sym_options(self, maxw, maxm):
        m, q = self.m, self.q
        out = []
        for fx in self.fixed_choices():
            for npairs in range(0, (maxw - len(fx)) // 2 + 1):
                for P in itertools.combinations(self.reps, npairs):
                    slots = list(fx) + list(P)
                    for mags in itertools.product((q, q - 2), repeat=len(slots)):
                        cm = sum(1 if s in self.fixed else 2 for s, mg in zip(slots, mags) if mg == q - 2)
                        if cm > maxm:
                            continue
                        for sg in itertools.product((1, -1), repeat=max(len(slots) - 1, 0)):
                            sg = (1,) + sg if slots else ()
                            Wd = {}
                            for s, mg, g in zip(slots, mags, sg):
                                Wd[s] = g * mg; Wd[(-s) % m] = g * mg
                            out.append((tuple(sorted(Wd.items())), len(Wd), cm))
        return out

    def free_options(self, maxw, maxm):
        m, q = self.m, self.q
        out = [((), 0, 0)]
        for j in range(1, maxw + 1):
            for P in itertools.combinations(range(1, m), j - 1):
                slots = (0,) + P
                for mags in itertools.product((q, q - 2), repeat=j):
                    cm = sum(1 for mg in mags if mg == q - 2)
                    if cm > maxm:
                        continue
                    for sg in itertools.product((1, -1), repeat=j - 1):
                        sg = (1,) + sg
                        out.append((tuple(sorted((s, g * mg) for s, mg, g in zip(slots, mags, sg))), j, cm))
        return out

    def norm_sym(self, Wt):
        if not Wt:
            return ()
        g = 1 if Wt[0][1] > 0 else -1
        return tuple((p, g * v) for p, v in Wt)

    def norm_free(self, Wt):
        if not Wt:
            return ()
        m = self.m
        best = None
        for p0, v0 in Wt:
            g = 1 if v0 > 0 else -1
            c = tuple(sorted(((p - p0) % m, g * v) for p, v in Wt))
            if best is None or c < best:
                best = c
        return best

    def canon(self, sk):
        m = self.m
        best = None
        for t in self.units:
            mt = [tuple(sorted(((t * p) % m, v) for p, v in Wt)) for Wt in sk]
            for sab in (0, 1):
                A, B = (mt[1], mt[0]) if sab else (mt[0], mt[1])
                A, B = self.norm_sym(A), self.norm_sym(B)
                for scd in (0, 1):
                    C, D = (mt[3], mt[2]) if scd else (mt[2], mt[3])
                    Cs = min(self.norm_free(C), self.norm_free(tuple(sorted(((-p) % m, v) for p, v in C))))
                    Ds = min(self.norm_free(D), self.norm_free(tuple(sorted(((-p) % m, v) for p, v in D))))
                    c = (A, B, Cs, Ds)
                    if best is None or c < best:
                        best = c
        return best

    def skeletons(self, nw, nm, raw=False):
        S = self.sym_options(nw, nm)
        F = self.free_options(nw, nm) if not raw else self.free_options_raw(nw, nm)
        if raw:
            S = [o for o in self.sym_options_raw(nw, nm)]
        byc = {}
        for o in F:
            byc.setdefault((o[1], o[2]), []).append(o[0])
        out = set()
        for a in S:
            for b in S:
                if a[1] + b[1] > nw or a[2] + b[2] > nm:
                    continue
                rw, rm = nw - a[1] - b[1], nm - a[2] - b[2]
                for (cw, cm), Cl in byc.items():
                    if cw > rw or cm > rm:
                        continue
                    Dl = byc.get((rw - cw, rm - cm), [])
                    for c in Cl:
                        for d in Dl:
                            sk = (a[0], b[0], c, d)
                            out.add(sk if raw else self.canon(sk))
        return sorted(out)

    def sym_options_raw(self, maxw, maxm):
        """All symmetric w-sets with all signs (no normalisation)."""
        out = []
        for Wt, nwk, cm in self.sym_options(maxw, maxm):
            out.append((Wt, nwk, cm))
            if Wt:
                out.append((tuple((p, -v) for p, v in Wt), nwk, cm))
        return out

    def free_options_raw(self, maxw, maxm):
        m, q = self.m, self.q
        out = [((), 0, 0)]
        for j in range(1, maxw + 1):
            for P in itertools.combinations(range(m), j):
                for mags in itertools.product((q, q - 2), repeat=j):
                    cm = sum(1 for mg in mags if mg == q - 2)
                    if cm > maxm:
                        continue
                    for sg in itertools.product((1, -1), repeat=j):
                        out.append((tuple(sorted((s, g * mg) for s, mg, g in zip(P, mags, sg))), j, cm))
        return out

    def random_skeleton(self, nw, nm, rng):
        """A random skeleton with N_w w-positions, N_- of them of magnitude q - 2: the block counts are drawn with
        probability proportional to the number of (position set, magnitude) choices, then positions, magnitudes and
        signs uniformly. Returns None if no skeleton has these counts."""
        m, q = self.m, self.q
        nrep = len(self.reps)
        fxs = self.fixed_choices()
        dists, weights = [], []
        for fA in fxs:
            for fB in fxs:
                for pA in range(0, nrep + 1):
                    for pB in range(0, nrep + 1):
                        rest = nw - len(fA) - len(fB) - 2 * pA - 2 * pB
                        if rest < 0:
                            continue
                        for nC in range(0, min(rest, m) + 1):
                            nD = rest - nC
                            if nD > m:
                                continue
                            slots = [1] * len(fA) + [2] * pA + [1] * len(fB) + [2] * pB + [1] * (nC + nD)
                            ways = [1] + [0] * nm          # magnitude choices with weighted minus count nm
                            for wgt in slots:
                                new = ways[:]
                                for c in range(nm + 1 - wgt):
                                    new[c + wgt] += ways[c]
                                ways = new
                            if ways[nm] == 0:
                                continue
                            cnt = (math.comb(nrep, pA) * math.comb(nrep, pB) * math.comb(m, nC) * math.comb(m, nD)
                                   * ways[nm])
                            dists.append((fA, fB, pA, pB, nC, nD)); weights.append(cnt)
        if not dists:
            return None
        w = np.array(weights, np.float64)
        fA, fB, pA, pB, nC, nD = dists[rng.choice(len(dists), p=w / w.sum())]
        reps = np.array(self.reps)
        while True:                                    # magnitudes: rejection on the minus count
            PA = sorted(rng.choice(reps, pA, replace=False).tolist()) if pA else []
            PB = sorted(rng.choice(reps, pB, replace=False).tolist()) if pB else []
            WC = sorted(rng.choice(m, nC, replace=False).tolist()) if nC else []
            WD = sorted(rng.choice(m, nD, replace=False).tolist()) if nD else []
            slotsA = list(fA) + PA
            slotsB = list(fB) + PB
            mags = {k: rng.choice([q, q - 2], len(sl)) for k, sl in enumerate((slotsA, slotsB, WC, WD))}
            cm = sum((1 if s_ in self.fixed else 2) * (mg == q - 2) for s_, mg in zip(slotsA, mags[0])) + \
                sum((1 if s_ in self.fixed else 2) * (mg == q - 2) for s_, mg in zip(slotsB, mags[1])) + \
                int(np.sum(mags[2] == q - 2)) + int(np.sum(mags[3] == q - 2))
            if cm == nm:
                break
        sk = []
        for k, sl in enumerate((slotsA, slotsB, WC, WD)):
            Wd = {}
            for s_, mg in zip(sl, mags[k]):
                v = int(mg) * int(rng.choice([1, -1]))
                Wd[s_] = v
                if k < 2:
                    Wd[(-s_) % m] = v
            sk.append(tuple(sorted(Wd.items())))
        return self.canon(tuple(sk))

    # ---------- stage 1 ----------
    def ylist(self, Wt, sym, raw=False):
        """All y on Z_m with the given w-positions and +-1 elsewhere (symmetric for A, B). Without raw, a block with
        no w-position has y(0) = +1 (negation)."""
        m = self.m
        Wd = dict(Wt)
        if sym:
            free = [a for a in self.reps if a not in Wd]
            ffree = [f for f in self.fixed if f not in Wd]       # 0 first when free
            S = allpm(len(free) + len(ffree))
            if not Wd and not raw:
                S = S[S[:, len(free)] == 1]                     # y(0) = +1 (negation)
            Y = np.zeros((len(S), m), np.int8)
            for p, v in Wd.items():
                Y[:, p] = v
            for j, a in enumerate(free):
                Y[:, a] = S[:, j]; Y[:, m - a] = S[:, j]
            for j, f in enumerate(ffree):
                Y[:, f] = S[:, len(free) + j]
        else:
            free = [a for a in range(m) if a not in Wd]
            S = allpm(len(free))
            if not Wd and not raw:
                S = S[S[:, 0] == 1]
            Y = np.zeros((len(S), m), np.int8)
            for p, v in Wd.items():
                Y[:, p] = v
            Y[:, free] = S
        return Y

    # ---------- canonical form of a stage-1 solution ----------
    def canon_y(self, ys):
        m = self.m
        rot = np.array([[(a + t) % m for a in range(m)] for t in range(m)])
        best = None
        for t in self.units:
            perm = np.array([(t * a) % m for a in range(m)])
            T = [y[perm] for y in ys]
            for sab in (0, 1):
                A, B = (T[1], T[0]) if sab else (T[0], T[1])
                A = A * (1 if A[0] > 0 else -1); B = B * (1 if B[0] > 0 else -1)
                for scd in (0, 1):
                    C, D = (T[3], T[2]) if scd else (T[2], T[3])
                    reps = []
                    for X in (C, D):
                        cands = np.concatenate([X[rot], X[::-1][rot]])
                        cands = cands * np.where(cands[:, :1] > 0, 1, -1)
                        # lexicographic minimum of the rows
                        o = np.lexsort(cands.T[::-1])
                        reps.append(tuple(cands[o[0]].tolist()))
                    c = (tuple(A.tolist()), tuple(B.tolist()), reps[0], reps[1])
                    if best is None or c < best:
                        best = c
        return best

    # ---------- stage 2 ----------
    def zlist(self, y, sym, fix_first=False):
        """All z for one block given its y: +-1 at the z-positions (|y| = 1), 0 at the w-positions; for A, B
        z(-a) = eps z(a)."""
        m, eps = self.m, self.eps
        zpos = [a for a in range(m) if abs(int(y[a])) == 1]
        if sym:
            free = [a for a in self.reps if abs(int(y[a])) == 1]
            ffree = [f for f in self.fixed if abs(int(y[f])) == 1] if eps == 1 else []
            S = allpm(len(free) + len(ffree))
            Z = np.zeros((len(S), m), np.int8)
            for j, a in enumerate(free):
                Z[:, a] = S[:, j]; Z[:, m - a] = eps * S[:, j]
            for j, f in enumerate(ffree):
                Z[:, f] = S[:, len(free) + j]
        else:
            S = allpm(len(zpos))
            if fix_first and len(zpos):
                S = S[S[:, 0] == 1]
            Z = np.zeros((len(S), m), np.int8)
            Z[:, zpos] = S
        return Z

    def uw(self, y):
        y = y.astype(np.int64)
        w = np.where(np.abs(y) > 1, np.sign(y), 0)
        u = y - self.q * w
        return u, w

    def stage2(self, ys, K=1, threads=None, normz=True):
        m, h, q, eps = self.m, self.h, self.q, self.eps
        U, Wv = zip(*[self.uw(y) for y in ys])
        R = np.zeros(h, np.int64)
        for u, w in zip(U, Wv):
            for s in range(1, h + 1):
                R[s - 1] += corr(u, w, s) + corr(w, u, s) + q * corr(w, w, s)
        # (III) is odd in s for eps = -1 (nothing to check at s = 0 or s = m/2), even for eps = +1
        S3 = [s_ for s_ in range(1, h + 1) if (2 * s_) % m] if eps == -1 else list(range(0, h + 1))
        # global z -> -z (a non-square multiplier of Z_q): fix the first free z of the largest list
        sizes = [int(np.sum(np.abs(y) == 1)) for y in ys]
        big = 2 + int(np.argmax(sizes[2:]))
        keys, Zs = [], []
        for k in range(4):
            Z = self.zlist(ys[k], k < 2, fix_first=(normz and k == big))
            Zi = Z.astype(np.int64)
            P = pafs(Zi, h)
            Cc = np.stack([eps * corr(U[k][None, :], Zi, s) + corr(Zi, U[k][None, :], s) for s in S3], axis=1)
            keys.append(np.concatenate([P, Cc], axis=1).astype(np.int16)); Zs.append(Z)
        tgt = np.concatenate([R + self.T2, np.zeros(len(S3), np.int64)]).astype(np.int16)
        total, sols = join(keys, tgt, K, threads)
        out = []
        for ids in sols:
            out.append([Zs[k][ids[k]] for k in range(4)])
        return total, out, U, Wv

    def lift(self, u, w, z):
        """The +-1 sequence of the block (u, w, z): on Z_N (q prime), or on Z_m x GF(q) in mixed-radix order."""
        if not self.prime:
            g = np.arange(self.N)
            a, yy = g // self.q, g % self.q
            x = np.where(yy == 0, u[a] + w[a], w[a] + z[a] * self.chi[yy])
            assert set(np.unique(x)) <= {-1, 1}
            return x.astype(np.int64)
        a, yy = self.ia, self.iy
        x = np.where(yy == 0, u[a] + w[a], w[a] + z[a] * self.chi[yy])
        assert set(np.unique(x)) <= {-1, 1}
        return x.astype(np.int64)


def record_hit(P, seqs, log):
    """Check a lifted solution as a complex Hadamard matrix and log it; returns the row sums."""
    rows = [int(x.sum()) for x in seqs]
    strs = ["".join('+' if v > 0 else '-' for v in x) for x in seqs]
    if P.bordered:
        import bsz
        Kr, Ki, b = bsz.build(seqs)
        rec = {"N": P.N, "m": P.v, "gens": [P.gen], "subgroup_order": (P.q - 1) // 2, "paley": [P.m, P.q],
               "rowsums": rows, "border": [[int(z.real), int(z.imag)] for z in b], "sequences": strs}
    elif P.prime:
        ssxx.complex_hadamard(seqs)
        rec = {"N": P.N, "gens": [P.gen], "subgroup_order": (P.q - 1) // 2, "paley": [P.m, P.q], "rowsums": rows,
               "sequences": strs}
    else:
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "verify"))
        import verify_complex as vc
        rec = {"N": P.N, "group": P.group, "field": {"q": P.q, "poly": P.poly},
               "generator": f"the squares of GF({P.q})*, trivial on Z_{P.m}", "subgroup_order": (P.q - 1) // 2,
               "paley": [P.m, P.q], "rowsums": rows, "sequences": strs}
        cert = {"order": 2 * P.N, "construction": "szollosi-ssxx-group", "group": P.group, "sequences": strs}
        vc.check(*vc.build(cert))
    if log:
        with open(log, 'a') as f:
            f.write(json.dumps(rec) + "\n")
    return rows


def join(keys, tgt, K, threads=None):
    sizes = [len(x) for x in keys]
    if min(sizes) == 0:
        return 0, []
    perm = min(([0, 1, 2, 3], [0, 2, 1, 3], [0, 3, 1, 2]),
               key=lambda o: sizes[o[0]] * sizes[o[1]] + sizes[o[2]] * sizes[o[3]])
    C = keys[0].shape[1]
    fd, path = tempfile.mkstemp(suffix='.bin', dir=W.tmpdir()); os.close(fd)
    env = dict(os.environ)
    if threads:
        env['OMP_NUM_THREADS'] = str(threads)
    try:
        with open(path, 'wb') as f:
            f.write(struct.pack('5i', *[sizes[k] for k in perm], C))
            f.write(np.asarray(tgt, np.int16).tobytes())
            for k in perm:
                f.write(np.ascontiguousarray(keys[k], dtype=np.int16).tobytes())
        out = subprocess.run([W.bin('join4x'), path, str(K)], capture_output=True, text=True, check=True,
                             env=env).stdout
    finally:
        os.remove(path)
    total = int(out.strip().split('TOTAL')[-1].split()[0])
    sols = []
    for line in out.splitlines():
        if line.startswith('SOL'):
            ids = list(map(int, line.split()[1:]))
            q = [None] * 4
            for slot, k in enumerate(perm):
                q[k] = ids[slot]
            sols.append(q)
    return total, sols


def run(N, q, all_=False, log=None, threads=None, maxsol=10 ** 6, raw=False, verbose=True, bordered=False):
    P = Paley(N, q, bordered)
    m, h = P.m, P.h
    t0 = time.time()
    adm = admissible(m, q, bordered)
    if P.eps == -1:
        adm = [(a, b) for a, b in adm if a >= 2 * len(P.fixed)]
    print(f"N={N}{' bordered' if bordered else ''}: {m} x {q}, eps={P.eps}, H=<{P.gen}> of order {(q - 1) // 2}; "
          f"admissible (N_w, N_-): {adm}",
          flush=True)
    found, nsk, n1, n1o, n2 = [], 0, 0, 0, 0
    for nw, nm in adm:
        sks = P.skeletons(nw, nm, raw=raw)
        print(f"  N_w={nw} N_-={nm}: {len(sks)} skeletons [{time.time() - t0:.0f}s]", flush=True)
        for sk in sks:
            nsk += 1
            lists = [P.ylist(sk[k], k < 2, raw=raw) for k in range(4)]
            keys = [pafs(Y, h).astype(np.int16) for Y in lists]
            tot, sols = join(keys, np.full(h, P.T1, np.int16), maxsol, threads)
            n1 += tot
            if tot > len(sols):
                print(f"  WARNING: {tot} stage-1 solutions, only {len(sols)} kept (--maxsol)", flush=True)
            if not sols:
                continue
            reps = set()                       # canonical stage-1 solutions as bytes (int8, 4 m entries)
            for ids in sols:
                ys = [lists[k][ids[k]] for k in range(4)]
                c = ys if raw else P.canon_y(ys)
                reps.add(np.asarray(c, np.int8).tobytes())
            n1o += len(reps)
            if verbose:
                print(f"  skeleton {sk}: stage 1 {tot} solutions, {len(reps)} orbits [{time.time() - t0:.0f}s]",
                      flush=True)
            for key in sorted(reps):
                ys = list(np.frombuffer(key, np.int8).astype(np.int64).reshape(4, m))
                tot2, zs, U, Wv = P.stage2(ys, K=1, threads=threads, normz=not raw)
                n2 += tot2
                if tot2 and (not found or (all_ and len(found) < 3)):
                    z = zs[0]
                    seqs = [P.lift(U[k], Wv[k], z[k].astype(np.int64)) for k in range(4)]
                    rows = record_hit(P, seqs, log)
                    print(f"  HIT N={N}: verified complex Hadamard matrix of order {2 * N}; rows {rows}; "
                          f"stage 2 count {tot2}", flush=True)
                    found.append(seqs)
                    if not all_ and not raw:
                        print(f"RESULT N={N} q={q}: FOUND [{time.time() - t0:.0f}s]", flush=True)
                        return found
    status = "FOUND" if found else "no solution (exhaustive)"
    print(f"RESULT N={N} q={q}: {status}; skeletons {nsk}, stage-1 solutions {n1} ({n1o} orbits), "
          f"stage-2 solutions {n2} [{time.time() - t0:.0f}s]", flush=True)
    return found


def sample(N, q, K, seed=1, log=None, threads=None, maxsol=10 ** 6, bordered=False):
    """Random skeletons (not exhaustive): stages 1 and 2 on each, stop at the first verified hit."""
    P = Paley(N, q, bordered)
    m, h = P.m, P.h
    rng = np.random.default_rng(seed)
    adm = [(a, b) for a, b in admissible(m, q, bordered) if P.eps == 1 or a >= 2 * len(P.fixed)]
    print(f"N={N}{' bordered' if bordered else ''}: {m} x {q}, eps={P.eps}: sampling {K} skeletons over {adm} "
          f"(seed {seed})", flush=True)
    t0, seen, n1 = time.time(), set(), 0
    for it in range(K):
        nw, nm = adm[rng.integers(len(adm))]
        sk = P.random_skeleton(nw, nm, rng)
        if sk is None or sk in seen:
            continue
        seen.add(sk)
        lists = [P.ylist(sk[k], k < 2) for k in range(4)]
        keys = [pafs(Y, h).astype(np.int16) for Y in lists]
        tot, sols = join(keys, np.full(h, P.T1, np.int16), maxsol, threads)
        if not sols:
            continue
        reps = set()
        for ids in sols:
            reps.add(np.asarray(P.canon_y([lists[k][ids[k]] for k in range(4)]), np.int8).tobytes())
        n1 += len(reps)
        print(f"  skeleton {it} ({nw},{nm}): stage 1 {tot} solutions, {len(reps)} orbits [{time.time() - t0:.0f}s]",
              flush=True)
        for key in reps:
            ys = list(np.frombuffer(key, np.int8).astype(np.int64).reshape(4, m))
            tot2, zs, U, Wv = P.stage2(ys, K=1, threads=threads)
            if tot2:
                seqs = [P.lift(U[k], Wv[k], zs[0][k].astype(np.int64)) for k in range(4)]
                rows = record_hit(P, seqs, log)
                print(f"  HIT N={N}: verified complex Hadamard matrix of order {2 * N}; rows {rows}; skeleton {sk}",
                      flush=True)
                print(f"RESULT N={N} q={q}: FOUND after {len(seen)} skeletons [{time.time() - t0:.0f}s]", flush=True)
                return True
    print(f"RESULT N={N} q={q}: nothing in {len(seen)} sampled skeletons ({n1} stage-1 orbits) [{time.time() - t0:.0f}s]",
          flush=True)
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("N", type=int); ap.add_argument("q", type=int)
    ap.add_argument("--all", action="store_true"); ap.add_argument("--raw", action="store_true")
    ap.add_argument("--log"); ap.add_argument("--threads", type=int)
    ap.add_argument("--maxsol", type=int, default=10 ** 6)
    ap.add_argument("--sample", type=int, help="random skeletons instead of all of them (stops at the first hit)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--bordered", action="store_true", help="group Z_m x GF(q) of order N - 1 (search/bsz.py)")
    a = ap.parse_args()
    if a.sample:
        sample(a.N, a.q, a.sample, a.seed, a.log, a.threads, a.maxsol, a.bordered)
    else:
        run(a.N, a.q, a.all, a.log, a.threads, a.maxsol, a.raw, bordered=a.bordered)


if __name__ == "__main__":
    main()

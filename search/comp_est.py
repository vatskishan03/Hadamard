#!/usr/bin/env python3
"""Compression-conditioned estimates for Szöllősi structures (search/ssxx.py, search/lowe_sz.py).

The Gaussian estimate of ssxx.gauss treats the class sums as independent continuous variables. Compressions show
that this is wrong in a structured way: for a divisor d of N, the fibre sums x(a) = sum_{g = a (mod d)} X(g) satisfy

    sum_X PAF_{x_X}(s) = 4N delta_s          on Z_d,

and when H acts on Z_d with large stabilisers the fibre sums lie on coarse progressions (for N = 3p and H trivial on
Z_3, x(a) = +-1 (mod 2|H|)), so most row-sum patterns have no solution at all and the others have far more than
the estimate. Here every block is split by its compression: the fibre classes are the H-orbits on Z_d (for the
symmetric blocks the <H, -1>-orbits), each block orbit lies over one fibre class and meets each of its fibres in
the same number of points, and the block's compression is a vector of fibre-class values. For a row-sum pattern:

  1. per block, every achievable compression and the number of H-invariant blocks with it (a product over fibre
     classes of one-dimensional counts);
  2. the compression tuples that satisfy the compressed equation exactly (a join on Z_d);
  3. per tuple, a multivariate lattice Gaussian for the class sums: blocks sampled uniformly with their
     compressions, mean and covariance of the class values, the lattice spanned by the samples (Hermite normal
     form), and the density at the target on the affine span.

E(pattern) = sum over feasible tuples of prod_X count_X(x_X) * P(class sums = target | compressions).

    python3 search/comp_est.py N g [g2 ...] --d 3 [--rows a,b,c,d] [--K 3000]
"""
import argparse, itertools, math, sys
from math import gcd, log10, pi
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ssxx
from bnw import Struct

KINDS = ssxx.KINDS                      # A, B symmetric ("L"), C, D arbitrary ("A")


class Comp:
    """fibre classes of Z_d under H (kind 'A') or <H, -1> (kind 'L'), and each block orbit's class and weight"""

    def __init__(self, st, d):
        assert st.m % d == 0
        self.st, self.d = st, d
        m = st.m
        Hd = sorted({h % d for h in st.H})
        Ld = sorted({h % d for h in st.L})
        self.fc, self.fid, self.w, self.orbs = {}, {}, {}, {}
        for kind, G in (("A", Hd), ("L", Ld)):
            cls, seen = [], {}
            for a in range(d):
                if a in seen:
                    continue
                O = sorted({a * g % d for g in G})
                for b in O:
                    seen[b] = len(cls)
                cls.append(O)
            self.fc[kind] = cls
            orb = st.orb(kind); sz = st.sizes(kind)
            fid, w = [], []
            for o in range(len(sz)):
                pts = np.flatnonzero(orb == o)
                f = seen[int(pts[0]) % d]
                assert all(seen[int(p) % d] == f for p in pts)
                w_ = len(pts) // len(cls[f])
                assert w_ * len(cls[f]) == len(pts)
                fid.append(f); w.append(w_)
            self.fid[kind], self.w[kind] = np.array(fid), np.array(w)

    def class_counts(self, kind):
        """per fibre class: dict value -> number of sign vectors on the orbits over it"""
        out = []
        for f in range(len(self.fc[kind])):
            cnt = {0: 1}
            for o in np.flatnonzero(self.fid[kind] == f):
                wo = int(self.w[kind][o]); nxt = {}
                for t, c in cnt.items():
                    nxt[t + wo] = nxt.get(t + wo, 0) + c
                    nxt[t - wo] = nxt.get(t - wo, 0) + c
                cnt = nxt
            out.append(cnt)
        return out

    def compressions(self, kind, r):
        """all achievable compressions (tuple of fibre-class values) with row sum r, with their counts"""
        cc = self.class_counts(kind)
        sizes = [len(c) for c in self.fc[kind]]
        out = {}
        for vals in itertools.product(*[sorted(c) for c in cc]):
            if sum(s * v for s, v in zip(sizes, vals)) != r:
                continue
            n = 1
            for c, v in zip(cc, vals):
                n *= c[v]
            out[vals] = n
        return out

    def expand(self, kind, vals):
        x = np.zeros(self.d, np.int64)
        for f, O in enumerate(self.fc[kind]):
            x[O] = vals[f]
        return x

    def sample(self, kind, vals, K, rng):
        """K uniform sign vectors (over the block orbits) with the given compression (vectorised backward DP)"""
        sz = self.st.sizes(kind); O = len(sz)
        Y = np.empty((K, O))
        for f, v in enumerate(vals):
            idx = np.flatnonzero(self.fid[kind] == f)
            ws = [int(self.w[kind][o]) for o in idx]
            n = len(ws); W_ = sum(ws)
            T = np.zeros((n + 1, 2 * W_ + 1))                  # T[i][t + W_]: ways to reach t with orbits i..
            T[n, W_] = 1.0
            for i in range(n - 1, -1, -1):
                w = ws[i]
                T[i, w:] += T[i + 1, :2 * W_ + 1 - w]            # +w
                T[i, :2 * W_ + 1 - w] += T[i + 1, w:]            # -w
            rem = np.full(K, v, np.int64)
            for i in range(n):
                w = ws[i]
                a = rem - w + W_; b = rem + w + W_
                cp = np.where((a >= 0) & (a <= 2 * W_), T[i + 1, np.clip(a, 0, 2 * W_)], 0.0)
                cm = np.where((b >= 0) & (b <= 2 * W_), T[i + 1, np.clip(b, 0, 2 * W_)], 0.0)
                plus = rng.random(K) * (cp + cm) < cp
                Y[:, idx[i]] = np.where(plus, 1.0, -1.0)
                rem = rem - np.where(plus, w, -w)
        return Y


def values(st, kind, Y):
    """class values y^T S_c y for all samples and classes at once"""
    S = st.S[kind].astype(np.float64)                         # (C, O, O)
    T = np.einsum("ko,cop->kcp", Y, S, optimize=True)
    return np.einsum("kcp,kp->kc", T, Y, optimize=True)


def cpaf(x):
    d = len(x)
    return np.array([int(np.dot(x, np.roll(x, -s))) for s in range(d)])


def comp_target(st, d):
    """the compressed target on Z_d: sum over the fibres of the class-sum target (0 plain, -4 bordered) off zero,
    and 4m at zero"""
    m, t = st.m, st.target
    tgt = np.full(d, t * (m // d), np.int64)
    tgt[0] = 4 * m + t * (m // d - 1)
    return tgt


def feasible_tuples(cp, N, rows):
    """compression tuples (one per block) with sum_X PAF_{x_X} = the compressed target on Z_d"""
    lists = [cp.compressions(k, r) for k, r in zip(KINDS, rows)]
    d = cp.d
    tgt = comp_target(cp.st, d)
    P = [{v: cpaf(cp.expand(k, v)) for v in L} for k, L in zip(KINDS, lists)]
    ab = {}
    for a in lists[0]:
        for b in lists[1]:
            key = tuple(P[0][a] + P[1][b])
            ab.setdefault(key, []).append((a, b))
    out = []
    for c in lists[2]:
        for dv in lists[3]:
            key = tuple(tgt - P[2][c] - P[3][dv])
            for a, b in ab.get(key, []):
                out.append((a, b, c, dv))
    return out, lists


def echelon(V):
    """integer echelon basis (rows) of the lattice spanned by the integer rows of V (unimodular row operations)"""
    rows = [list(map(int, r)) for r in V if any(r)]
    if not rows:
        return []
    ncol = len(rows[0])
    B = []
    piv = 0
    M = rows
    for c in range(ncol):
        # gcd-reduce column c among rows of M not yet used as pivots
        cand = [r for r in M if r[c] != 0]
        rest = [r for r in M if r[c] == 0]
        while len(cand) > 1:
            cand.sort(key=lambda r: abs(r[c]))
            p = cand[0]
            new = [p]
            for r in cand[1:]:
                q = r[c] // p[c]
                r2 = [x - q * y for x, y in zip(r, p)]
                if r2[c] != 0:
                    new.append(r2)
                elif any(r2):
                    rest.append(r2)
            cand = new
        if cand:
            B.append(cand[0])
        M = rest
    return B


def bareiss_det(M):
    """exact determinant of a square integer matrix (fraction-free elimination)"""
    M = [list(r) for r in M]
    n = len(M); sign = 1; prev = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if M[i][k] != 0), None)
            if sw is None:
                return 0
            M[k], M[sw] = M[sw], M[k]; sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]


def hnf_covolume(V):
    """covolume of the lattice spanned by the integer rows of V (in its real span) and its rank; the Gram
    determinant is computed exactly (the echelon basis can have very large entries)"""
    B = echelon(V)
    if not B:
        return 1.0, 0, B
    G = [[sum(x * y for x, y in zip(r1, r2)) for r2 in B] for r1 in B]
    det = bareiss_det(G)
    return math.exp(0.5 * math.log(det)), len(B), B


def in_lattice(B, v):
    """is the integer vector v in the lattice with echelon basis B?"""
    v = list(map(int, v))
    for r in B:
        c = next(i for i, x in enumerate(r) if x)
        if v[c] % r[c]:
            return False
        q = v[c] // r[c]
        v = [x - q * y for x, y in zip(v, r)]
    return not any(v)


def block_stats(st, cp, kind, vals, K, rng):
    """cached per (kind, compression): count, mean and covariance of the class values, one value vector (the
    lattice base point) and an echelon basis of the lattice spanned by the differences"""
    cache = cp.__dict__.setdefault("stats_cache", {})
    key = (kind, vals)
    if key not in cache:
        Y = cp.sample(kind, vals, K, rng)
        V = values(st, kind, Y)
        Vi = np.rint(V).astype(np.int64)
        U = np.unique(Vi, axis=0)
        cache[key] = (cp.compressions_cache[key], V.mean(0), np.cov(V.T), Vi[0], echelon((U[1:] - U[0])[:120]))
    return cache[key]


def cond_estimate(st, cp, tup, K, rng, target=None):
    """log10 of P(class sums = target | compressions) under the multivariate lattice Gaussian, the log10 of the
    number of block quadruples with these compressions, the rank of the covariance and of the lattice"""
    tv = st.target if target is None else target
    t = np.full(st.C, float(tv))
    mu = np.zeros(st.C); S = np.zeros((st.C, st.C)); base = np.zeros(st.C, np.int64); gens = []
    ln = 0.0
    for kind, vals in zip(KINDS, tup):
        cnt, m_, c_, b_, B_ = block_stats(st, cp, kind, vals, K, rng)
        ln += log10(cnt); mu += m_; S += c_; base += b_; gens += B_
    B = echelon(gens) if gens else []
    if B:
        G = [[sum(x * y for x, y in zip(r1, r2)) for r2 in B] for r1 in B]
        covol, rank = math.exp(0.5 * math.log(bareiss_det(G))), len(B)
    else:
        covol, rank = 1.0, 0
    if float(tv) == int(tv):
        if not in_lattice(B, np.full(st.C, int(tv), np.int64) - base):
            return -math.inf, ln, -1, rank
    w, Uv = np.linalg.eigh(S)
    keep = w > 1e-7 * max(w.max(), 1e-12)
    delta = t - mu
    proj = Uv[:, keep] @ (Uv[:, keep].T @ delta)
    if np.linalg.norm(delta - proj) > 1e-6 * (1 + np.linalg.norm(delta)):
        return -math.inf, ln, int(keep.sum()), rank
    k_ = int(keep.sum())
    q = float(np.sum((Uv[:, keep].T @ delta) ** 2 / w[keep]))
    lp = log10(covol) - 0.5 * k_ * log10(2 * pi) - 0.5 * float(np.sum(np.log10(w[keep]))) - q / 2 / math.log(10)
    return lp, ln, k_, rank


def estimate(N, gens, d, rows=None, K=3000, seed=1, verbose=True, plain=True):
    st = Struct(N, gens, plain=plain)
    cp = Comp(st, d)
    cp.compressions_cache = {}
    rng = np.random.default_rng(seed)
    if rows:
        pats = [rows]
    elif plain:
        pats = ssxx.patterns(st)
    else:
        import bsz
        pats = [r for _, r in bsz.patterns(st)]
    res = []
    for pat in pats:
        tups, lists = feasible_tuples(cp, N, pat)
        for kind, L in zip(KINDS, lists):
            for v, n in L.items():
                cp.compressions_cache[(kind, v)] = n
        tot = 0.0; best = None
        for tup in tups:
            lp, ln, k_, rank = cond_estimate(st, cp, tup, K, rng)
            if lp == -math.inf:
                continue
            e = lp + ln
            tot += 10 ** e
            if best is None or e > best[0]:
                best = (e, tup, ln, lp, k_, rank)
        if plain:
            g0 = ssxx.gauss(st, pat)
        else:
            import bsz
            g0 = bsz.gauss(st, pat)
        le = log10(tot) if tot > 0 else -math.inf
        res.append((le, pat, len(tups), best, g0))
        if verbose:
            print(f"N={N} H={st.gens} d={d} rows {tuple(pat)}: {len(tups)} feasible compression tuples, "
                  f"conditioned logE {le:.2f} (plain Gaussian {g0:.2f})"
                  + (f"; best tuple {best[1]} logE {best[0]:.2f} (pools 10^{best[2]:.2f}, P 10^{best[3]:.2f}, "
                     f"rank {best[4]}/{best[5]})" if best else ""), flush=True)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("N", type=int); ap.add_argument("gens", type=int, nargs="+")
    ap.add_argument("--d", type=int, required=True); ap.add_argument("--rows")
    ap.add_argument("--K", type=int, default=3000); ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--bordered", action="store_true")
    a = ap.parse_args()
    rows = tuple(int(x) for x in a.rows.split(",")) if a.rows else None
    estimate(a.N, a.gens, a.d, rows, a.K, a.seed, plain=not a.bordered)


if __name__ == "__main__":
    main()

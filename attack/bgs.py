"""Bordered Goethals-Seidel arrays from cyclotomic products.

**The bordered array.** Let G be an abelian group, K a subgroup of index q with Q = G/K, and T a set of t characters
of Q that is a union of Galois orbits. Put Gamma = 4 * sum_{chi in T} chi, an integer function on Q with
Gamma(0) = 4t. If four +-1 sequences on G satisfy

    sum_i PAF_i(g) = -Gamma(g + K)        for every g != 0,

and the border problem below is solvable, then the Goethals-Seidel array M developed over G, bordered by 4t rows and
columns that are constant on the cosets of K, is a Hadamard matrix of order 4n with n = |G| + t:

    H = [ X  R ]        R[x, (i, g)] = Rh[x, i, g + K],   C[(i, g), x] = Ch[i, g + K, x].
        [ C  M ]

The border problem only involves the compressions y_i of the sequences to Q (their sums over the cosets of K):
X is a Hadamard matrix of order 4t, Rh (4t x 4q) has orthogonal rows, Ch Ch^T = I_4 (x) Gamma, and
X Ch^T + Rh Mc^T = 0, where Mc is the Goethals-Seidel array of the y_i over Q. (Spence's construction is Q = 1,
T = {1}; the matrices of Alpoge et al. for n = 3 (mod 4) use Q = Z_2 x Z_2 with T the three nontrivial characters.)
The fibres of the compressed problem satisfy sum_i |y_i^(chi)|^2 = 4t for chi in T and 4n for chi not in T.

**Cyclotomic products.** Here G = Z_m x GF(r) with r a prime not dividing m, K = qZ_m x GF(r) (q | m), and T the
characters of Z_q whose orders lie in a set D of divisors of q (t = sum_{d in D} phi(d), Gamma = 4 sum c_d with
Ramanujan sums c_d). The sequences are invariant under H = {h : h mod m in K_m, h mod r in C_0} where C_0 is the
subgroup of index e of GF(r)^* and K_m = <kgen> <= Z_m^*. A block is given at each a in Z_m by its value f(a) at
(a, 0) and its values x_0(a) .. x_{e-1}(a) on the cyclotomic classes (a, C_j); its compression to Z_m is
z(a) = f(a) + F s(a), F = (r - 1)/e, s(a) = sum_j x_j(a). The compressed targets on Z_m are
T_1(0) = 4mr - 4t(r - 1) and T_1(b) = -r Gamma(b mod q), and summing squares gives

    sum over the 4m positions of  F s^2 + 2 f s  =  4 (m - t) e,

so for t close to m only a few spikes (positions with s != 0) fit, and the large negative targets -r Gamma(b) are
met by pairs of spikes of opposite sign at distance b.

  python3 attack/bgs.py scan N [N ...] [--tmax T] [--maxlist L]     structures with stage-1 solutions for 4N
"""
import argparse, itertools, json, math, os, struct, subprocess, sys, tempfile, time
from math import comb, gcd
import numpy as np
from pathlib import Path as _P
sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
import verify


def primitive_root(p):
    t = p - 1; fs = [f for f in range(2, t + 1) if t % f == 0 and all(f % d for d in range(2, int(f ** 0.5) + 1))]
    return next(r for r in range(2, p) if all(pow(r, t // f, p) != 1 for f in fs))


def is_prime(p):
    return p > 1 and all(p % d for d in range(2, int(p ** 0.5) + 1))


def phi(n):
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)


def mobius(n):
    res, p = 1, 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            if n % p == 0: return 0
            res = -res
        p += 1
    return -res if n > 1 else res


def ramanujan(d, a):
    g = gcd(a, d); return mobius(d // g) * phi(d) // phi(d // g)


def gamma_vec(q, D):
    """Gamma(a) = 4 sum_{d in D} c_d(a) on Z_q."""
    return np.array([4 * sum(ramanujan(d, a % q) for d in D) for a in range(q)], dtype=np.int64)


class Structure:
    def __init__(self, m, r, e, q, D, kgen=1):
        assert is_prime(r) and m % r and (r - 1) % e == 0 and e % 2 == 0 and m % q == 0
        assert all(q % d == 0 for d in D)
        self.m, self.r, self.e, self.q, self.D = m, r, e, q, tuple(sorted(D))
        self.t = sum(phi(d) for d in D)
        self.n = m * r + self.t
        self.F = (r - 1) // e
        self.Gamma = gamma_vec(q, D)
        assert self.Gamma[0] == 4 * self.t
        g = primitive_root(r); self.g = g
        ind = np.zeros(r, dtype=np.int64); x = 1
        for k in range(r - 1): ind[x] = k; x = x * g % r
        self.cls = np.full(r, -1, dtype=np.int64); self.cls[1:] = ind[1:] % e
        self.half = int(self.cls[r - 1])                    # class of -1
        self.nb = e // 2 if self.half else e
        self.reps = [0] + [int(np.flatnonzero(self.cls == c)[0]) for c in range(e)]
        self.kgen = kgen % m if m > 1 else 1
        K = {1}
        while True:
            K2 = K | {k * self.kgen % m for k in K}
            if K2 == K: break
            K = K2
        self.K = sorted(K)
        self.korb = np.full(m, -1, dtype=np.int64); self.kreps = []
        for a in range(m):
            if self.korb[a] < 0:
                for k in self.K: self.korb[a * k % m] = len(self.kreps)
                self.kreps.append(a)
        self.T1 = np.array([4 * m * r - 4 * self.t * (r - 1)] + [-r * int(self.Gamma[b % q]) for b in range(1, m)])
        self.budget = 4 * (m - self.t) * e
        v = m * r; self.v = v
        ge = pow(g, e, r)
        self.gens = [next(h for h in range(1, v) if h % m == 1 and h % r == ge)]
        if len(self.K) > 1: self.gens.append(next(h for h in range(1, v) if h % m == self.kgen and h % r == 1))

    def describe(self):
        return (f"n={self.n} = {self.m}*{self.r} + {self.t}: q={self.q} D={list(self.D)} e={self.e} F={self.F} "
                f"K_m=<{self.kgen}> |H|={self.F * len(self.K)}")

    def target(self, gidx):
        """Target of sum PAF at the element of Z_v with index gidx (CRT: a = g mod m)."""
        if gidx == 0: return 4 * self.v
        return -int(self.Gamma[(gidx % self.m) % self.q])

    # pattern ids as in cyclo_product: id = (fbit << e) | mask
    def pattern(self, pid):
        e = self.e
        return 1 - 2 * ((pid >> e) & 1), [1 - 2 * ((pid >> j) & 1) for j in range(e)]

    def patterns(self, f, s):
        e = self.e; fb = 0 if f == 1 else 1
        return [(fb << e) | mk for mk in range(1 << e) if e - 2 * bin(mk).count('1') == s]

    def sequence(self, slices):
        v, m, r = self.v, self.m, self.r; t = np.arange(v)
        P = np.array([[self.pattern(p)[0]] + self.pattern(p)[1] for p in slices])
        col = np.where(t % r == 0, 0, 1 + self.cls[t % r])
        return P[t % m, col]


def stage1(S, maxpairs=4e6, maxopts=4e5, maxspikes=6):
    """All quadruples of compressions z_i on Z_m (K_m-invariant, z = f + F s) with sum_i PAF_{z_i} = T_1, one per
    class under negation, translation (K_m trivial) and reversal of single blocks."""
    m, e, F = S.m, S.e, S.F
    B = S.budget
    if B < 0: return []
    cost = lambda s, f: F * s * s + 2 * f * s
    spikes = [(f, s) for s in range(-e, e + 1, 2) if s for f in (1, -1) if cost(s, f) <= B]
    orbits = [[a for a in range(m) if S.korb[a] == o] for o in range(len(S.kreps))]
    units = [len(o) for o in orbits]
    if len(orbits) > 22: raise RuntimeError(f'stage 1: {len(orbits)} position orbits')
    if spikes and B // min(cost(s, f) for f, s in spikes) > maxspikes:
        raise RuntimeError(f'stage 1: up to {B // min(cost(s, f) for f, s in spikes)} spikes')
    out = [((), ())]
    def rec(start, used, cur, c):
        for i in range(start, len(units)):
            for f, s in spikes:
                cc = c + units[i] * cost(s, f)
                if cc > B: continue
                out.append((used + (i,), cur + ((f, s),)))
                if len(out) > maxopts: raise RuntimeError('stage 1: too many spike placements')
                rec(i + 1, used + (i,), cur + ((f, s),), cc)
    rec(0, (), (), 0)
    Y, C = [], []
    for pos, fsl in out:
        rest = [i for i in range(len(orbits)) if i not in pos]
        c0 = sum(units[i] * cost(s, f) for i, (f, s) in zip(pos, fsl))
        if len(Y) > maxopts: raise RuntimeError('stage 1: too many options')
        for signs in itertools.product((1, -1), repeat=len(rest)):
            y = [0] * m
            assign = dict(zip(pos, fsl)); assign.update({i: (sg, 0) for i, sg in zip(rest, signs)})
            for i, (f, s) in assign.items():
                for a in orbits[i]: y[a] = f + F * s
            Y.append(y); C.append(c0)
    Y = np.array(Y, dtype=np.int64).reshape(-1, m); C = np.array(C, dtype=np.int64)
    P = np.stack([(Y * np.roll(Y, -b, axis=1)).sum(1) for b in range(1, m // 2 + 1)], axis=1) if m > 1 else np.zeros((len(Y), 0), np.int64)
    tgt = S.T1[1:m // 2 + 1]
    def canon_ids(Y):
        best = None
        for sgn in (1, -1):
            for rv in (1, -1):
                for t in (range(m) if len(S.K) == 1 else [0]):
                    Z = sgn * Y[:, (rv * np.arange(m) + t) % m]
                    if best is None: best = Z.copy(); continue
                    lt = np.zeros(len(Y), bool); eq = np.ones(len(Y), bool)
                    for a in range(m):
                        lt |= eq & (Z[:, a] < best[:, a]); eq &= Z[:, a] == best[:, a]
                    best[lt] = Z[lt]
        u, inv = np.unique(best, axis=0, return_inverse=True)
        return inv.ravel()
    cid = canon_ids(Y)
    R = np.random.default_rng(5).integers(1, 2 ** 61, size=max(1, m // 2))
    hk = (P * R[:P.shape[1]]).sum(1); htgt = int((tgt * R[:len(tgt)]).sum()) if len(tgt) else 0
    buckets = {int(c): np.flatnonzero(C == c) for c in np.unique(C)}
    found = {}
    for c0, i0 in buckets.items():
        for c1, i1 in buckets.items():
            if c1 < c0: continue
            for c2, i2 in buckets.items():
                if c2 < c1: continue
                c3 = B - c0 - c1 - c2
                if c3 < c2 or c3 not in buckets: continue
                i3 = buckets[c3]
                if len(i0) * len(i1) > maxpairs or len(i2) * len(i3) > maxpairs:
                    raise RuntimeError(f'stage 1: {len(i0) * len(i1)} x {len(i2) * len(i3)} pairs')
                A_, B_ = np.meshgrid(i0, i1, indexing='ij'); A_, B_ = A_.ravel(), B_.ravel()
                k01 = hk[A_] + hk[B_]; o = np.argsort(k01); ks = k01[o]
                Cg, Dg = np.meshgrid(i2, i3, indexing='ij'); Cg, Dg = Cg.ravel(), Dg.ravel()
                need = htgt - (hk[Cg] + hk[Dg])
                lo = np.searchsorted(ks, need, 'left'); hi = np.searchsorted(ks, need, 'right')
                for rr in np.flatnonzero(hi > lo):
                    for p_ in range(lo[rr], hi[rr]):
                        a, b, c, d = A_[o[p_]], B_[o[p_]], Cg[rr], Dg[rr]
                        if np.any(P[a] + P[b] + P[c] + P[d] - tgt): continue
                        key = tuple(sorted((cid[a], cid[b], cid[c], cid[d])))
                        if key not in found:
                            found[key] = tuple(tuple(int(x) for x in Y[j]) for j in (a, b, c, d))
    return [found[k] for k in sorted(found)]


def list_size(S, y):
    """Members of a block with compression y (product over K_m-orbits of the class patterns)."""
    tot = 1
    for a in S.kreps:
        f, s = fs_of(S, y[a])
        tot *= comb(S.e, (S.e + s) // 2)
    return tot


def fs_of(S, yv):
    s = int(round(yv / S.F)); f = yv - S.F * s
    assert f in (1, -1), (yv, S.F)
    return f, s


def nequations(S):
    """Number of shift classes (b, C_c) with c a class (not y = 0), under K_m, C_0 and -1."""
    m, e = S.m, S.e
    seen = set(); cnt = 0
    for b in range(m):
        for c in range(e):
            if (b, c) in seen: continue
            cnt += 1
            todo = [(b, c)]
            while todo:
                x = todo.pop()
                if x in seen: continue
                seen.add(x)
                bb, cc = x
                for k in S.K: todo.append((bb * k % m, cc))
                todo.append(((-bb) % m, (cc + S.half) % e))
    return cnt


def candidate_structures(n, tmax=40, rmin=11, mmax=30):
    """(m, r, t, q, D) with m r = n - t, r prime, q | m, q >= t, D a set of divisors of q with sum phi = t."""
    out = []
    for t in range(1, tmax + 1):
        v = n - t
        for r in range(rmin, v + 1):
            if v % r or not is_prime(r): continue
            m = v // r
            if m % r == 0 or m < t or m > mmax: continue
            for q in range(t, m + 1):
                if m % q: continue
                divs = [d for d in range(1, q + 1) if q % d == 0]
                for k in range(1, len(divs) + 1):
                    for Dset in itertools.combinations(divs, k):
                        # rows of Rh and Ch are +-1 functions with Fourier support in T: for odd q the sum over a
                        # block is odd, so the trivial character must lie in T
                        if q % 2 and 1 not in Dset: continue
                        if sum(phi(d) for d in Dset) == t:
                            out.append((m, r, t, q, Dset))
    return out


def scan(ns, tmax=40, maxlist=3e7, log=print, maxspikes=4):
    for n in ns:
        cands = candidate_structures(n, tmax)
        log(f"n={n}: {len(cands)} (m, r, t, q, D) candidates")
        for m, r, t, q, D in cands:
            for e in range(2, r, 2):
                if (r - 1) % e: continue
                F = (r - 1) // e
                if F < 2: continue
                # quick necessary condition: a spike costs at least F*4 - 4 per position
                budget = 4 * (m - t) * e
                for kgen in sorted({k for k in range(1, m) if gcd(k, m) == 1} or {1}):
                    try:
                        S = Structure(m, r, e, q, D, kgen)
                    except AssertionError:
                        continue
                    if len(S.K) > 1 and any(S.K == Structure(m, r, e, q, D, k2).K for k2 in range(1, kgen) if gcd(k2, m) == 1):
                        continue
                    # filters: few spikes, and lists that can be small
                    cmin = 4 * F - 4
                    if budget and budget // cmin > maxspikes: continue
                    if budget and budget < cmin: continue
                    no = len(S.kreps)
                    if comb(e, e // 2) ** max(0, no - budget // cmin) > 1e3 * maxlist: continue
                    try:
                        quads = stage1(S)
                    except RuntimeError as ex:
                        log(f"  {S.describe()}: {ex}"); continue
                    if not quads: continue
                    sizes = [[list_size(S, y) for y in qd] for qd in quads]
                    best = min(max(sz) for sz in sizes)
                    log(f"  {S.describe()}: {len(quads)} stage-1 classes, eqs {nequations(S)}, "
                        f"smallest max list {best:.3g}" + ("" if best <= maxlist else " (too big)"))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    sc = sub.add_parser('scan'); sc.add_argument('n', type=int, nargs='+'); sc.add_argument('--tmax', type=int, default=40)
    sc.add_argument('--maxlist', type=float, default=3e7)
    a = ap.parse_args()
    log = lambda s: print(s, flush=True)
    if a.cmd == 'scan': scan(a.n, a.tmax, a.maxlist, log)

"""Gaussian estimates for bordered Goethals-Seidel structures (NOTES section 18).

The sequences are constant on the orbits of a group S of automorphisms of G (spence_lib.Structure). For a row-sum
pattern (r_1, .., r_4) with sum r_i^2 = 4|G| + sum_{d != 0} target(d), each block is sampled uniformly with its
row sum, its class values (autocorrelations on the classes of <S, -1>) are reduced to their mean, variance and lattice,
and the four blocks are combined as independent Gaussians on each class (dropping the class the row sums determine).
log10 E = sum log10(pool sizes) + sum over classes of the log-density at the target. The model ignores the exact
compression constraints and can be wrong by many orders of magnitude in either direction (NOTES section 18).

  python3 search/bgs_gauss.py z2z3 N K1,K2,..    family Q = Z_2^2 x Z_3, T = (nontrivial on Z_2^2) x (all on Z_3):
                                                 G = Z_2^2 x Z_3 x Z_K1 x ..., one coupled generator
                                                 (alpha, multipliers on K) per choice of multiplier subgroups
"""
import itertools, math, sys, time
from math import gcd, log10, pi
from pathlib import Path as _P
import numpy as np
sys.path.insert(0, str(_P(__file__).resolve().parent))
import spence_lib as L


def row_patterns(total, parity):
    out = []; lim = math.isqrt(total)
    rs = [r for r in range(0, lim + 1) if (r - parity) % 2 == 0]
    for a in rs:
        for b in rs:
            if b < a: continue
            for c in rs:
                if c < b: continue
                d2 = total - a * a - b * b - c * c
                if d2 < c * c: continue
                d = math.isqrt(d2)
                if d * d == d2 and (d - parity) % 2 == 0: out.append((a, b, c, d))
    return out


def estimate(st, gamma, samples=1200, seed=1):
    """(sum of squares of the row sums, [(log10 E, pattern)] best first, log10 of the total)."""
    v = st.G.v
    tgt = -np.array([gamma[st.reps[c]] for c in range(st.C)], dtype=float)
    tot = 4 * v - int(sum(int(gamma[d]) for d in range(1, v)))
    rng = np.random.default_rng(seed); cache = {}
    def stats(r):
        if r not in cache:
            N = st.pool(r)
            if N == 0: cache[r] = None
            else:
                V = st.class_values(st.sample(r, samples, rng)); Vi = np.rint(V).astype(np.int64); lat = []
                for c in range(st.C):
                    u = np.unique(Vi[:, c]); g = int(np.gcd.reduce(np.diff(u))) if len(u) > 1 else 0
                    lat.append((g, int(u[0] % g) if g else int(u[0])))
                cache[r] = (N, V.mean(0), V.var(0), lat)
        return cache[r]
    res = []
    for pat in row_patterns(tot, v % 2):
        S4 = [stats(r) for r in pat]
        if any(s is None for s in S4): continue
        mu = sum(s[1] for s in S4); var = sum(s[2] for s in S4); terms = []; ok = True
        for c in range(st.C):
            g = 0; resid = 0
            for s in S4: gg, aa = s[3][c]; g = gcd(g, gg); resid += aa
            if var[c] < 1e-9:
                if abs(tgt[c] - mu[c]) > 1e-6: ok = False; break
                continue
            if g and (tgt[c] - resid) % g: ok = False; break
            terms.append((st.csize[c] * math.sqrt(var[c]),
                          -0.5 * log10(2 * pi * var[c]) - (tgt[c] - mu[c]) ** 2 / (2 * var[c]) / math.log(10) + log10(max(g, 1))))
        if not ok: continue
        terms.sort()
        res.append((sum(log10(s[0]) for s in S4) + sum(t for _, t in terms[:-1]), pat))
    res.sort(reverse=True)
    total = (res[0][0] + log10(sum(10 ** (e - res[0][0]) for e, _ in res))) if res else -math.inf
    return tot, res, total


def unit_reps(m):
    if m <= 2: return [1]
    seen = {}
    for a in (a for a in range(1, m) if gcd(a, m) == 1):
        H = frozenset(pow(a, j, m) for j in range(1, m + 1))
        seen.setdefault(H, a)
    return list(seen.values())


def z2z3(n, kmods, maxorb=90, minorb=30):
    mods = [2, 2, 3] + kmods; G = L.AbGroup(mods); co = G.coords
    assert G.v == n - 9, "need n = 9 + 12 |K|"
    u0 = (co[:, 0] == 0) & (co[:, 1] == 0)
    gamma = np.where(co[:, 2] == 0, np.where(u0, 36, -12), 0)
    AL = np.array([[0, 1], [1, 1]], dtype=np.int64)
    for combo in itertools.product(*[unit_reps(m) for m in kmods]):
        A = np.eye(len(mods), dtype=np.int64); A[0:2, 0:2] = AL
        for i, mu in enumerate(combo): A[3 + i, 3 + i] = mu
        st = L.Structure(G, [G.perm_from_matrix(A)])
        if not minorb <= st.O <= maxorb: continue
        t = time.time(); tot, res, total = estimate(st, gamma)
        print(f"n={n} K={kmods} sigma=(alpha, {combo}): orbits {st.O}, classes {st.C}, sum r^2 = {tot}, "
              f"log10 E = {total:.1f} (best pattern {res[0][1] if res else None}) [{time.time() - t:.0f}s]", flush=True)


if __name__ == '__main__':
    if sys.argv[1] == 'z2z3':
        z2z3(int(sys.argv[2]), [int(x) for x in sys.argv[3].split(',')])

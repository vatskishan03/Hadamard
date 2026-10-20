"""Goethals-Seidel structures over arbitrary abelian groups G (|G| = v odd), symmetry H <= Aut(G).

Four +-1 sequences constant on H-orbits with sum_i P_i(g) = 0 for g != 0 give a Hadamard matrix of
order 4v via the Goethals-Seidel array (group-developed blocks, R: g -> -g).  The row sums satisfy
sum r_i^2 = 4v.  Same Gaussian model as search/triage.py (per-class lattice, one class dropped)."""
import math
from math import gcd, log10, pi
import numpy as np


def patterns(v):
    out, lim = [], int(math.isqrt(4 * v)) + 1
    odd = range(1, lim + 1, 2)
    for a in odd:
        for b in odd:
            if b < a: continue
            for c in odd:
                if c < b: continue
                d2 = 4 * v - a * a - b * b - c * c
                if d2 < c * c: continue
                d = math.isqrt(d2)
                if d * d == d2: out.append((a, b, c, d))
    return out


def gs_estimate(st, samples=3000, log_list=8.0, seed=1, max_patterns=None):
    v = st.G.v
    assert v % 2 == 1
    pm = v % 4
    target = -pm
    pats = [q for q in patterns(v) if all(st.pool(r) > 0 for r in q)]
    if not pats:
        return None
    rng = np.random.default_rng(seed)
    stats = {}
    for r in sorted({x for q in pats for x in q}):
        V = (st.class_values(st.sample(r, samples, rng)) - pm) / 4.0
        Vi = np.rint(V).astype(np.int64)
        lat = []
        for c in range(st.C):
            u = np.unique(Vi[:, c])
            g = int(np.gcd.reduce(np.diff(u))) if len(u) > 1 else 0
            lat.append((g, int(u[0] % g) if g else int(u[0])))
        stats[r] = (V.mean(0), V.var(0), lat)
    best = None
    for q in pats:
        mu = sum(stats[r][0] for r in q); var = sum(stats[r][1] for r in q)
        terms, dead = [], False
        for c in range(st.C):
            g, res = 0, 0
            for r in q:
                gg, aa = stats[r][2][c]; g, res = gcd(g, gg), res + aa
            if var[c] < 1e-9:
                if abs(target - mu[c]) > 1e-6: dead = True; break
                continue
            if g and (target - res) % g: dead = True; break
            t = -0.5 * log10(2 * pi * var[c]) - (target - mu[c]) ** 2 / (2 * var[c]) / math.log(10) + log10(max(g, 1))
            terms.append((st.csize[c] * math.sqrt(var[c]), t))
        if dead: continue
        terms.sort()
        lE = sum(log10(st.pool(r)) for r in q) + sum(t for _, t in terms[:-1])
        lN = log10(min(st.pool(r) for r in q))
        if lN <= log_list: lW = 2 * lN - lE
        else: lW = 3 * lN - lE - log_list - 0.5 * log10(6 / pi) * (st.C - 1)
        if best is None or lW < best["lW"]:
            best = dict(pattern=q, lE=lE, lN=lN, lW=lW)
    return best

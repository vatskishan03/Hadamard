"""Triage skew Goethals-Seidel structures: A skew-type (A(0)=1, A(-x)=-A(x)) constant on H-orbits (-1 not in H),
B, C, D constant on H-orbits; sum of PAF = 4n delta gives a skew-Hadamard matrix of order 4n (Wallis-Whiteman).
Gaussian model for E (per-class lattice, one class dropped) and the k-tree work W ~ N_big^3 / (E L) style estimate."""
import sys, math
from math import gcd, log10, pi
import numpy as np
sys.path.insert(0, '/home/user/Hadamard/attack')
import spence_lib as L
from gs_sweep_subgroups import subgroups

def skew_sampler(st, n, count, rng):
    # orbits paired by negation: pair index for each orbit
    neg = st.G.neg()
    partner = np.array([st.orb[neg[st.reps_orb[o]]] for o in range(st.O)])
    return partner

def analyse(n, gens, samples=2500, rng=None):
    G = L.AbGroup([n]); st = L.Structure(G, [G.perm_from_matrix(np.array([[g]])) for g in gens])
    rep = np.zeros(st.O, dtype=np.int64)
    for o in range(st.O): rep[o] = np.flatnonzero(st.orb == o)[0]
    partner = st.orb[G.neg()[rep]]
    if np.any(partner[1:] == np.arange(1, st.O)): return None      # an orbit fixed by -1: no skew A
    pairs = [(o, partner[o]) for o in range(1, st.O) if o < partner[o]]
    NA = 2 ** len(pairs)
    # skew A samples
    YA = np.zeros((samples, st.O)); YA[:, 0] = 1
    s = np.where(rng.random((samples, len(pairs))) < 0.5, 1.0, -1.0)
    for j, (o, q) in enumerate(pairs): YA[:, o] = s[:, j]; YA[:, q] = -s[:, j]
    # admissible row sums for B, C, D: 1 + rB^2 + rC^2 + rD^2 = 4n
    pats = []
    lim = int(math.isqrt(4 * n))
    for b in range(1, lim + 1, 2):
        for c in range(b, lim + 1, 2):
            d2 = 4 * n - 1 - b * b - c * c
            if d2 < c * c: continue
            d = math.isqrt(d2)
            if d * d == d2 and all(st.pool(r) > 0 for r in (b, c, d)): pats.append((b, c, d))
    if not pats: return None
    pm = n % 4
    VA = st.class_values(YA)
    best = None
    for q in pats:
        Vs = [VA] + [st.class_values(st.sample(r, samples, rng)) for r in q]
        V = sum(Vs)
        lp = []; dead = False
        for c in range(st.C):
            us = [np.unique(v[:, c]) for v in Vs]
            g = 0
            for u in us:
                if len(u) > 1: g = gcd(g, int(np.gcd.reduce(np.diff(u).astype(np.int64))))
            res = int(sum(u[0] for u in us))
            if g and (0 - res) % g: dead = True; break
            var = V[:, c].var(); mu = V[:, c].mean()
            if var < 1e-9:
                if abs(mu) > 1e-6: dead = True; break
                continue
            lp.append((st.csize[c] * math.sqrt(var), -0.5 * log10(2 * pi * var) - mu ** 2 / (2 * var) / math.log(10) + log10(max(g, 1))))
        if dead: continue
        lp.sort()
        lE = log10(NA) + sum(log10(st.pool(r)) for r in q) + sum(t for _, t in lp[:-1])
        lNmax = max(log10(st.pool(r)) for r in q)
        if best is None or lE > best[0]: best = (lE, q, log10(NA), lNmax)
    return (st.O, st.C, best) if best else None

if __name__ == '__main__':
    ns = [int(a) for a in sys.argv[2:]]; maxo = int(sys.argv[1])
    rng = np.random.default_rng(1)
    for n in ns:
        for H, gens in subgroups(n).items():
            if (n - 1) in H: continue
            seen = set(); O = 0
            for t in range(n):
                if t not in seen: seen |= {t * h % n for h in H}; O += 1
            if O > maxo or O < 6: continue
            r = analyse(n, list(gens), rng=rng)
            if r is None: continue
            O, C, (lE, q, lNA, lNmax) = r
            print(f"n={n} 4n={4*n} |H|={len(H)} gens={gens} orbits={O} classes={C} rows=(1,{q}) log10E={lE:.2f} log10NA={lNA:.1f} log10Nmax={lNmax:.1f}", flush=True)

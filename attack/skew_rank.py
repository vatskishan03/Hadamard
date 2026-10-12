"""Rank skew Goethals-Seidel multiplier structures by the yield of a k-tree search on sub-sampled lists.

For each n, each subgroup H of Z_n^* with -1 not in H, and each row-sum pattern (1; b, c, d), the Gaussian
model of skew_triage.py gives the expected number E of solutions. Lists capped at L distinct sequences
keep a fraction prod_i min(1, L / N_i) of them:  E' = E * prod_i min(1, L / N_i)  (A included).

  python3 skew_rank.py L MAXO n ...
"""
import math, sys
from math import gcd, log10, pi
import numpy as np
sys.path.insert(0, '/home/user/Hadamard/attack')
import spence_lib as L
from gs_sweep_subgroups import subgroups
GAIN = 0.21   # log10 gain per class from low-energy B, C, D (spread 2.4x narrower), measured on v = 205


def patterns(st, n):
    lim = int(math.isqrt(4 * n)); out = []
    for b in range(1, lim + 1, 2):
        for c in range(b, lim + 1, 2):
            d2 = 4 * n - 1 - b * b - c * c
            if d2 < c * c: continue
            d = math.isqrt(d2)
            if d * d == d2 and all(st.pool(r) > 0 for r in (b, c, d)): out.append((b, c, d))
    return out


def rank(n, gens, cap, samples=2500, rng=None):
    G = L.AbGroup([n]); st = L.Structure(G, [G.perm_from_matrix(np.array([[g]])) for g in gens])
    return st, rank_structure(st, n, cap, samples, rng)


def rank_structure(st, n, cap, samples=2500, rng=None):
    G = st.G
    rep = np.array([np.flatnonzero(st.orb == o)[0] for o in range(st.O)])
    partner = st.orb[G.neg()[rep]]
    if np.any(partner[1:] == np.arange(1, st.O)): return []
    pairs = [(o, partner[o]) for o in range(1, st.O) if o < partner[o]]
    lNA = len(pairs) * log10(2)
    YA = np.zeros((samples, st.O)); YA[:, 0] = 1
    sg = np.where(rng.random((samples, len(pairs))) < 0.5, 1.0, -1.0)
    for j, (o, q) in enumerate(pairs): YA[:, o] = sg[:, j]; YA[:, q] = -sg[:, j]
    VA = st.class_values(YA)
    cache = {}
    def cv(r):
        if r not in cache: cache[r] = st.class_values(st.sample(r, samples, rng))
        return cache[r]
    out = []
    for q in patterns(st, n):
        Vs = [VA] + [cv(r) for r in q]
        V = sum(Vs); lp = []; dead = False
        for c in range(st.C):
            us = [np.unique(v[:, c]) for v in Vs]
            g = 0
            for u in us:
                if len(u) > 1: g = gcd(g, int(np.gcd.reduce(np.rint(np.diff(u)).astype(np.int64))))
            res = int(round(sum(u[0] for u in us)))
            if g and (0 - res) % g: dead = True; break
            var = V[:, c].var(); mu = V[:, c].mean()
            if var < 1e-9:
                if abs(mu) > 1e-6: dead = True; break
                continue
            lp.append((st.csize[c] * math.sqrt(var), -0.5 * log10(2 * pi * var) - mu ** 2 / (2 * var) / math.log(10) + log10(max(g, 1))))
        if dead: continue
        lp.sort()
        lN = [lNA] + [log10(st.pool(r)) for r in q]
        lE = sum(lN) + sum(t for _, t in lp[:-1])
        lE2 = lE + sum(min(0.0, log10(cap) - x) for x in lN)
        out.append((lE2, lE, q, lN))
    return out


if __name__ == '__main__':
    cap = float(sys.argv[1]); maxo = int(sys.argv[2]); ns = [int(a) for a in sys.argv[3:]]
    rng = np.random.default_rng(1)
    for n in ns:
        for H, gens in subgroups(n).items():
            if (n - 1) in H: continue
            O = len({min(t * h % n for h in H) for t in range(n)})
            if O > maxo or O < 6: continue
            st, res = rank(n, list(gens), cap, rng=rng)
            for lE2, lE, q, lN in sorted(res, reverse=True)[:2]:
                lG = lE2 + GAIN * (st.C - 1)        # low-energy lists (skew_gen) for B, C, D; A uniform
                if max(lE2, lG) < -3: continue
                print(f"{lE2:7.2f} gen {lG:6.2f} n={n} 4n={4*n} |H|={len(H)} gens={','.join(map(str, gens))} orbits={st.O} classes={st.C} "
                      f"rows=(1,{q[0]},{q[1]},{q[2]}) log10E={lE:.2f} log10N=({','.join(f'{x:.1f}' for x in lN)})", flush=True)

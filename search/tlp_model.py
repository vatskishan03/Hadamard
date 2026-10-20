#!/usr/bin/env python3
"""Twisted Legendre pairs over Z_l: (A, A.mu) with A constant on H-orbits, mu^2 in <H,-1>, mu not in <H,-1>.
Condition: P_A(s) + P_A(mu s) = -2 on the classes of <H, mu, -1>.  Gaussian estimate of E (one class
dropped, per-class lattice), pool N (row sum 1).  usage: tlp_model.py MAXO n [n...]"""
import os
import sys, math
from math import gcd, log10, pi
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lp_search as LP

def orbit_labels(l, H):
    lab = -np.ones(l, dtype=np.int64); k = 0
    for t in range(l):
        if lab[t] < 0:
            for h in H: lab[t * h % l] = k
            k += 1
    return lab, k

def twisted_structures(l, maxo):
    subs = LP.subgroups(l)
    units = [u for u in range(1, l) if gcd(u, l) == 1]
    out = []
    for H, gens in subs.items():
        lab, k = orbit_labels(l, H)
        if k > maxo: continue
        L = frozenset(H) | frozenset((-h) % l for h in H)
        seen = set()
        for mu in units:
            if mu in L or (mu * mu % l) not in L: continue
            coset = frozenset(mu * x % l for x in L)
            if coset in seen: continue
            seen.add(coset)
            out.append((H, gens, mu, lab, k))
    return out

def model(l, H, mu, lab, k, samples=4000, rng=None):
    sizes = np.bincount(lab, minlength=k)
    K = set(H) | {mu * h % l for h in H}
    KL = K | {(-x) % l for x in K}
    cls = -np.ones(l, dtype=np.int64); reps = []
    for s in range(1, l):
        if cls[s] < 0:
            for x in KL: cls[s * x % l] = len(reps)
            reps.append(s)
    C = len(reps)
    # pool size with row sum 1 (DP over orbit sizes)
    cnt = {0: 1}
    for s in sizes:
        nxt = {}
        for t, c in cnt.items():
            nxt[t + s] = nxt.get(t + s, 0) + c; nxt[t - s] = nxt.get(t - s, 0) + c
        cnt = nxt
    N = cnt.get(1, 0)
    if N == 0: return None
    # sample uniformly among row-sum-1 sequences by rejection-free DP sampling (approx: rejection)
    rng = rng or np.random.default_rng(1)
    Y = []
    while sum(len(y) for y in Y) < samples:
        Z = np.where(rng.random((20000, k)) < 0.5, 1, -1)
        Y.append(Z[(Z @ sizes) == 1])
    Y = np.concatenate(Y)[:samples]
    x = np.arange(l)
    A = Y[:, lab]                                   # (S, l)
    V = np.zeros((len(Y), C))
    for c, s in enumerate(reps):
        V[:, c] = (A * A[:, (x + s) % l]).sum(1) + (A * A[:, (x + mu * s) % l]).sum(1)
    terms = []
    for c in range(C):
        u = np.unique(V[:, c]); g = int(np.gcd.reduce(np.diff(u).astype(np.int64))) if len(u) > 1 else 0
        mu_c, var = V[:, c].mean(), V[:, c].var()
        if var < 1e-9:
            if abs(mu_c + 2) > 1e-6: return dict(lE=-math.inf, N=N, C=C, k=k)
            continue
        if g and (int(-2 - u[0])) % g: return dict(lE=-math.inf, N=N, C=C, k=k)
        terms.append((sizes.max() * math.sqrt(var), -0.5 * log10(2 * pi * var) - (mu_c + 2) ** 2 / (2 * var) / math.log(10) + log10(max(g, 1))))
    terms.sort()
    lE = log10(N) + sum(t for _, t in terms[:-1])
    return dict(lE=lE, N=N, C=C, k=k)

if __name__ == '__main__':
    maxo = int(sys.argv[1]); ns = [int(a) for a in sys.argv[2:]]
    rng = np.random.default_rng(1)
    for n in ns:
        l = 2 * n - 1
        rows = []
        for H, gens, mu, lab, k in twisted_structures(l, maxo):
            if k < 8: continue
            m = model(l, H, mu, lab, k, rng=rng)
            if m is None: continue
            rows.append((m['lE'], len(H), gens, mu, k, m['C'], m['N']))
        rows.sort(key=lambda r: -r[0] if math.isfinite(r[0]) else 1e9)
        tot = sum(10 ** r[0] for r in rows if math.isfinite(r[0]))
        print(f"n={n} l={l}: {len(rows)} twisted structures (orbits <= {maxo}); sum E = {tot:.3g}", flush=True)
        for r in rows[:6]:
            print(f"   lE={r[0]:7.2f} |H|={r[1]} gens={r[2]} mu={r[3]} orbits={r[4]} classes={r[5]} pool={r[6]:.3g}", flush=True)

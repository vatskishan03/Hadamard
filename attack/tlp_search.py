#!/usr/bin/env python3
"""Exhaustive search for twisted Legendre pairs (A, A.mu) over Z_l, l = 2n - 1.
For every multiplier subgroup H with at most MAXO orbits and every mu (mu^2 in <H,-1>, mu not in it, one per
coset), enumerate all A constant on H-orbits with row sum 1 and test P_A(s) + P_A(mu s) = -2 on all classes of
<H, mu, -1>.  A hit gives a Legendre pair (A, A.mu) of length l and H(2l + 2) = H(4n).
usage: tlp_search.py MAXO MAXPOOL n [n ...]      (n = 'all' for every open order)"""
import sys, re, json, time, math
import numpy as np
sys.path.insert(0, '/home/user/Hadamard/attack')
import lp_search as LP, tlp_model as TM

def class_matrix(l, lab, k, reps, mu):
    S = np.zeros((len(reps), k, k))
    x = np.arange(l)
    for c, s in enumerate(reps):
        np.add.at(S[c], (lab, lab[(x + s) % l]), 1)
        np.add.at(S[c], (lab, lab[(x + mu * s) % l]), 1)
    return S

def search(l, H, mu, lab, k, maxpool):
    sizes = np.bincount(lab, minlength=k)
    K = set(H) | {mu * h % l for h in H}
    KL = K | {(-x) % l for x in K}
    cls = -np.ones(l, dtype=np.int64); reps = []
    for s in range(1, l):
        if cls[s] < 0:
            for x in KL: cls[s * x % l] = len(reps)
            reps.append(s)
    S = class_matrix(l, lab, k, reps, mu)
    hits = []; total = 0
    if 2 ** k > maxpool * 64: return None, 0
    for start in range(0, 1 << k, 1 << 20):
        idx = np.arange(start, min(start + (1 << 20), 1 << k), dtype=np.int64)
        Y = 1 - 2 * ((idx[:, None] >> np.arange(k)) & 1)
        Y = Y[(Y @ sizes) == 1].astype(np.float64)
        if len(Y) == 0: continue
        total += len(Y)
        V = np.einsum('bo,cop,bp->bc', Y, S, Y, optimize=True)
        ok = np.all(np.abs(V + 2) < 0.5, axis=1)
        for i in np.flatnonzero(ok): hits.append(Y[i].astype(int))
    return hits, total

def verify_pair(l, a, b):
    for s in range(1, l):
        if int(a @ np.roll(a, -s) + b @ np.roll(b, -s)) != -2: return False
    return abs(int(a.sum())) == 1 and abs(int(b.sum())) == 1

if __name__ == '__main__':
    maxo, maxpool = int(sys.argv[1]), float(sys.argv[2])
    if sys.argv[3] == 'all':
        ns = [int(m.group(2)) for line in open('/home/user/Hadamard/FRONTIER.md') for m in [re.match(r'\|\s*(\d+)\s*\|\s*(\d+)\s*\|', line)] if m]
    else:
        ns = [int(a) for a in sys.argv[3:]]
    for n in ns:
        l = 2 * n - 1; t0 = time.time(); tried = 0; covered = 0
        for H, gens, mu, lab, k in TM.twisted_structures(l, maxo):
            hits, total = search(l, H, mu, lab, k, maxpool)
            if hits is None: continue
            tried += 1; covered += total
            for y in hits:
                a = y[lab]; b = a[(np.arange(l) * pow(mu, -1, l)) % l]   # b(x) = a(mu^-1 x): b = A.mu^{-1}; any image works
                b2 = a[(np.arange(l) * mu) % l]
                ok = verify_pair(l, a, b2) or verify_pair(l, a, b)
                print(f"*** n={n} l={l} |H|={len(H)} gens={gens} mu={mu} orbits={k}: twisted LP found, verified={ok}", flush=True)
                with open('/home/user/Hadamard/attack/tlp_hits.json', 'a') as f:
                    f.write(json.dumps(dict(n=n, l=l, H=sorted(H), gens=list(gens), mu=mu,
                                            a=''.join('+' if v > 0 else '-' for v in a), verified=bool(ok))) + "\n")
        print(f"n={n} l={l}: {tried} twisted structures searched, {covered:.3g} sequences ({time.time()-t0:.0f}s)", flush=True)

"""Instance-level completion estimator (HEURISTIC): mu = prod_X |F_X| * (Gaussian local-limit density of the sum of four independent uniform fibre keys at 0, lattice spacing per component),
with the fibres F_X built EXACTLY (all 2^nv sign vectors of the free t signs, tau-set and PSD filters as in gl2.c) and keys = (pafS + q PAF_t at the I lags, Q at the Q lags).
mu_of_instance(sh, S, E) -> (mu, fibre sizes) is the library function.
Driver (K = Z43, q = 11, H_K = <6>, the slot (0,4,7,10)/(0,4,6,10) of the hit of skew/1892, the level-1 sample of the hit round): the instances are given to level 2 (skf), and mu is computed for
the completing instances and a random sample of 25 others.  Observed: 16 completing instances of the 4032 against an expected ~0.03 by mu (mean mu ~ 7e-6 per instance, about 600 times too low):
the completions are far from random key coincidences; no available estimator predicts which K-symmetric shape is rich (HEURISTIC).
usage (from search/skewplug):  python3 estimates/inst_mu.py        (a few minutes)"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
DATA = os.path.join(HERE, '..', 'data'); RUNS = os.path.join(HERE, '..', 'runs')
import json, math, random, time
from sk import SShape
from sk_pat import row_patterns_skew
from sk_samp import level1
from gl1 import L1
from skf_io import skf_header, run_skf
import numpy as np, itertools
from math import gcd


def free_orbits(sh, b, e):
    if b == 0:
        fl = sh.A_dflags(e)
        return [sh.at_orbits[o] for o in range(len(sh.at_orbits)) if not fl[o]]
    K = sh.kCD
    return [[(int(p), int(K.sig[p])) for p in K.pos[i]] for i in range(K.nu) if e[int(K.pos[i][0])] == 0]


def fibre(sh, b, s, e, nI_lags, Qlags, W, T):
    orbs = free_orbits(sh, b, e); nv = len(orbs); N = sh.N
    M = 1 << nv
    bits = ((np.arange(M)[:, None] >> np.arange(nv)[None, :]) & 1)
    sg = (1 - 2 * bits).astype(np.int8)                   # (M, nv)
    t = np.zeros((M, N), dtype=np.int8)
    for j, orb in enumerate(orbs):
        for (p, sig) in orb: t[:, p] = sg[:, j] * sig
    s = s.astype(np.int64)
    # PSD filter: |sum_z (s + i sqrt(q) t) w^(psi z)|^2 <= T + 0.01, all psi = 1..N-1 (W: (N-1, N) complex)
    sc = W @ s.astype(np.complex128)                      # (N-1,)
    tc = t.astype(np.float64) @ W.T                       # (M, N-1)
    val = sc[None, :] + 1j * math.sqrt(sh.r) * tc
    ok = (np.abs(val) ** 2 <= T + 0.01).all(1)
    return t, ok, nv


def keys(sh, t, s, lagsI, lagsQ):
    N = sh.N; q = sh.r; eps = sh.eps; s = s.astype(np.int64)
    tt = t.astype(np.int32)
    cols = []
    for k in lagsI:
        ipp = sh.shift[k]
        pafS = int((s * s[ipp]).sum())
        cols.append(pafS + q * (tt * tt[:, ipp]).sum(1))
    for k in lagsQ:
        w = s[sh.shift[k]] + eps * s[sh.shneg[k]]
        cols.append(tt @ w.astype(np.int32))
    return np.stack(cols, 1).astype(np.int64)


def mu_of_instance(sh, S, E, W=None):
    N = sh.N; T = sh.T2
    if W is None: W = np.exp(2j * np.pi * np.outer(np.arange(1, N), np.arange(N)) / N)
    lagsI = [int(k) for k in sh.repsPM if k != 0]; lagsQ = [int(k) for k in sh.repsQ]
    Ks = []; sizes = []
    for b in range(4):
        t, ok, nv = fibre(sh, b, S[b], E[b], lagsI, lagsQ, W, T)
        t = t[ok]; sizes.append(len(t))
        if len(t) == 0: return 0.0, sizes
        Ks.append(keys(sh, t, S[b], lagsI, lagsQ))
    d = Ks[0].shape[1]
    mean = sum(K.mean(0) for K in Ks); cov = sum(np.cov(K.T) if len(K) > 1 else np.zeros((d, d)) for K in Ks)
    # lattice spacing per component and coset test
    g = np.zeros(d, dtype=np.int64); off = np.zeros(d, dtype=np.int64)
    for K in Ks:
        for j in range(d):
            dif = np.abs(K[:, j] - K[0, j]); gj = int(np.gcd.reduce(dif)) if len(dif) else 0
            g[j] = gcd(int(g[j]), gj); off[j] += K[0, j]
    for j in range(d):
        if g[j] == 0:
            if off[j] != 0: return 0.0, sizes
            g[j] = 1
        elif off[j] % g[j] != 0: return 0.0, sizes
    var = np.diag(cov); keep = var > 1e-9
    if (~keep).any():
        if np.any(np.abs(mean[~keep]) > 1e-6): return 0.0, sizes
    c = cov[np.ix_(keep, keep)] + 1e-6 * np.eye(keep.sum()); m = mean[keep]; sp = g[keep].astype(float)
    sign, logdet = np.linalg.slogdet(c); z = -m
    q = float(z @ np.linalg.solve(c, z))
    lg = np.log(sp).sum() - 0.5 * keep.sum() * math.log(2 * math.pi) - 0.5 * logdet - 0.5 * q
    logmu = lg + sum(math.log(max(s_, 1)) for s_ in sizes)
    return math.exp(logmu) if logmu < 700 else float('inf'), sizes


if __name__ == '__main__':
    sh = SShape([43], 11, [6]); l1 = L1(sh)
    al = (0, 4, 7, 10); ms = (0, 4, 6, 10); rnd = 9
    slots = []
    for l in open(os.path.join(DATA, 'cost473.jsonl')):
        r = json.loads(l)
        if 't_l2_inst' in r and r.get('n1_per_round', 0) > 0:
            c1 = r['t_l1_round'] / r['n1_per_round']; eff = 10 ** r['logrho'] / (c1 + r['t_l2_inst']); slots.append((math.log10(eff), tuple(r['al']), tuple(r['ms']), r.get('logN1', 0), r.get('logrho', 0)))
    slots.sort(key=lambda x: -x[0]); slots = slots[:40]; wts = [10 ** (s[0] - slots[0][0]) for s in slots]
    rng = random.Random(1 * 1000003 + rnd); rng.choices(slots, weights=wts, k=1); pats = row_patterns_skew(sh, al, ms); A = sh.enum_A(al[0], ms[0]); pat = rng.choice(pats)
    cnt, S, E, st = level1(sh, l1, al, ms, pat, 10000, 60000, 60000, seed=rng.randrange(10 ** 6), A_list=A, maxout=20000, hseed=rnd + 1)
    print('instances', cnt)
    res, ss = run_skf(sh, [[(S[i, b], E[i, b]) for b in range(4)] for i in range(len(S))], maxsol=1000, tlimit=60.0, hdr=skf_header(sh))
    comp = np.array([res[i]['compl'] for i in range(len(S))])
    idx = list(np.flatnonzero(comp > 0)[:3]) + random.Random(1).sample(range(len(S)), 25)
    W = np.exp(2j * np.pi * np.outer(np.arange(1, sh.N), np.arange(sh.N)) / sh.N)
    mus = []; t0 = time.time()
    for i in idx:
        mu, sizes = mu_of_instance(sh, S[i], E[i], W); mus.append(mu)
        print(i, 'observed completions', int(comp[i]), 'mu %.3g' % mu, 'fibres', sizes, 'skf fib', res[i]['fib'], flush=True)
    print('time per instance %.2f s' % ((time.time() - t0) / len(idx)))
    rs = [m for i, m in zip(idx, mus) if comp[i] == 0]
    print('mean mu over non-completing sample %.3g; expected completions in all %d instances ~ %.3g (observed total completions %d)' % (np.mean(rs), len(S), np.mean(rs) * len(S), comp.sum()))

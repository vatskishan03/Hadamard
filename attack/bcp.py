"""Bordered cyclotomic products (attack/bgs.py, NOTES §18-19) with the low-energy machinery of cp_lowe.py.

The structure is bgs.Structure: G = Z_m x GF(r), K = qZ_m x GF(r), T a set of characters of Z_q, targets
sum_i PAF_i(g) = -Gamma(g + K) for g != 0. The level-B/level-C keys of cp_lowe are the same; only their targets
change: level C at (b, d_c) is -Gamma(b mod q), level B at d_c is -(m/q) sum_alpha Gamma(alpha).

  python3 attack/bcp.py scan n1,n2,.. [--tmax 15] [--mmax 60] [--maxspikes 12] [--samples 60000] [--L 1e7]
      stage 1 (exact, capped) and Gaussian estimates for every candidate structure -> $HADAMARD_WORK/bcp_scan.jsonl
"""
import argparse, json, math, sys, time
import numpy as np
from math import comb, gcd
from pathlib import Path as _P
sys.path[:0] = [str(_P(__file__).resolve().parent)]
import workdir as W
import cyclo_product as CP
import cp_lowe as CL
import bgs


class Adapter:
    """bgs.Structure seen through the interface of cyclo_product.Structure (plain, field size r)."""

    def __init__(self, Sb):
        self.Sb = Sb
        self.v, self.q, self.e, self.m, self.F, self.skew = Sb.v, Sb.r, Sb.e, Sb.m, Sb.F, False
        self.cls, self.half, self.nb, self.reps, self.g = Sb.cls, Sb.half, Sb.nb, Sb.reps, Sb.g
        self.K, self.korb, self.kreps, self.kgen = Sb.K, Sb.korb, Sb.kreps, Sb.kgen
        self.h = Sb.gens[0]; self.bkey = []; self.dkey = Sb.nb

    pattern = CP.Structure.pattern
    patterns = CP.Structure.patterns
    cc_table = CP.Structure.cc_table
    sequence = CP.Structure.sequence


def targets(Sb, D):
    """Key targets: level B (nb coordinates) then level C ((m // 2) e coordinates)."""
    m, q, e = Sb.m, Sb.q, Sb.e
    G = [int(Sb.Gamma[a % q]) for a in range(m)]
    tB = [-sum(G)] * Sb.nb
    tC = [-G[b] for b in range(1, m // 2 + 1) for c in range(e)]
    t = np.array(tB + tC, dtype=float)
    assert len(t) == D
    return t


def estimate(Sb, y, samples, L, rng):
    """(log10 E of the class, log10 E inside low-energy lists of L members, rank, sizes)."""
    S = Adapter(Sb); O = CL.Orbits(S, 0)
    P = CL.prepare(S, O, y, samples, rng)
    t = targets(Sb, P['D'])
    lists = [b[2] for b in P['blocks']]
    est, rank, off = CL.gauss_estimate(lists, t)
    est += sum(math.log10(b[0].size / len(b[2])) for b in P['blocks'])
    # inside lists: members below the energy quantile, with the targets t
    mus, sigs = P['mus'], P['sigs']
    Sp = np.linalg.pinv(sum(sigs), rcond=1e-10); mu = sum(mus)
    kept, corr = [], 0.0
    for k in range(4):
        blk, I, K = P['blocks'][k]
        ck = mus[k] + sigs[k] @ Sp @ (t - mu); Wk, _ = CL.whitening(sigs[k])
        en = (((K - ck) @ Wk.T) ** 2).sum(1)
        frac = min(1.0, L / blk.size)
        sel = en <= (float(np.quantile(en, frac)) if frac < 1 else float(en.max() + 1))
        kept.append(K[sel]); corr += math.log10(min(L, blk.size) / sel.sum())
    estL = CL.gauss_estimate(kept, t)[0] + corr
    return est, estL, rank, off, P['sizes']


def scan(ns, tmax, mmax, maxspikes, samples, L, out):
    done = set()
    if out.exists():
        for l in open(out):
            r = json.loads(l); done.add((r['n'], r['m'], r['r'], r['t'], r['q'], tuple(r['D']), r['e'], r['kgen']))
    for n in ns:
        for m, r, t, q, D in bgs.candidate_structures(n, tmax, rmin=5, mmax=mmax):
            for e in range(2, r, 2):
                if (r - 1) % e or (r - 1) // e < 2: continue
                seenK = []
                for kgen in sorted({k for k in range(1, m) if gcd(k, m) == 1} or {1}):
                    try: Sb = bgs.Structure(m, r, e, q, D, kgen)
                    except AssertionError: continue
                    if Sb.K in seenK: continue
                    seenK.append(Sb.K)
                    key = (n, m, r, t, q, tuple(D), e, kgen)
                    if key in done: continue
                    F = Sb.F; cmin = 4 * F - 4
                    if Sb.budget < 0 or (Sb.budget and Sb.budget < cmin): continue
                    mem0 = comb(e, e // 2) ** len(Sb.kreps)
                    if mem0 > 1e13 or mem0 < 1e4: continue
                    base = dict(n=n, m=m, r=r, t=t, q=q, D=list(D), e=e, F=F, kgen=kgen, K=len(Sb.K), reps=len(Sb.kreps),
                                mem0=round(math.log10(mem0), 1))
                    t0 = time.time()
                    try:
                        quads = bgs.stage1(Sb, maxpairs=3e7, maxopts=2e6, maxspikes=maxspikes)
                    except (RuntimeError, MemoryError) as ex:
                        rec = dict(base, error=str(ex)[:80], sec=round(time.time() - t0))
                        print(json.dumps(rec), flush=True)
                        with open(out, 'a') as f: f.write(json.dumps(rec) + '\n')
                        continue
                    rec = dict(base, classes=len(quads))
                    if quads:
                        best = None
                        for ci, y in enumerate(quads[:4]):
                            try:
                                est, estL, rank, off, sizes = estimate(Sb, y, samples, L, np.random.default_rng(1))
                            except Exception as ex:
                                continue
                            if best is None or est > best[0]:
                                best = (est, estL, rank, off, sizes, ci, y)
                        if best:
                            est, estL, rank, off, sizes, ci, y = best
                            rec.update(config=ci, y=y, logE=round(est, 2), logE_lists=round(estL, 2), rank=rank, off=round(off, 2),
                                       sizes=[round(math.log10(max(1, x)), 1) for x in sizes])
                    rec['sec'] = round(time.time() - t0)
                    print(json.dumps(rec), flush=True)
                    with open(out, 'a') as f: f.write(json.dumps(rec) + '\n')


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    sc = sub.add_parser('scan'); sc.add_argument('ns')
    sc.add_argument('--tmax', type=int, default=15); sc.add_argument('--mmax', type=int, default=60)
    sc.add_argument('--maxspikes', type=int, default=12); sc.add_argument('--samples', type=int, default=60000)
    sc.add_argument('--L', type=float, default=1e7)
    a = ap.parse_args()
    if a.cmd == 'scan':
        scan([int(x) for x in a.ns.split(',')], a.tmax, a.mmax, a.maxspikes, a.samples, a.L, W.WORK / 'bcp_scan.jsonl')


if __name__ == '__main__':
    main()

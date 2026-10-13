"""Sweep of cyclotomic-product structures for the low-energy k-tree of cp_lowe.py.

For every target v (skew: odd v of the skew census; plain: n of FRONTIER.md), every factorisation v = m q with q prime
not dividing m and m >= 2, every even e | q - 1 (skew: (q - 1)/e odd), every cyclic K_m = <k> <= Z_m^* (skew: no
position orbit closed under negation) and every admissible twist t (the generator acts as (a, y) -> (k a, g^t y);
|K_m| t = 0 mod e; skew: t a multiple of 2^v2(e), so that the class of -1 stays outside the fibre symmetry), whose
non-skew blocks have between MINM and MAXM members, stage 1 is run (cached by cp_lowe.stage1_cached) and, for up to
NCLS stage-1 classes, the Gaussian estimate of the number of solutions is computed for the whole class and inside
low-energy lists of L members per non-skew block. One JSON line per class:

  python3 attack/cp_sweep.py skew|plain v1,v2,.. [--L 1e7] [--ncls 6] [--minm 1e5] [--maxm 3e10]

When the exact stage 1 is too large, NCLS classes are sampled with cp_lowe.stage1_sample (seed 1; 'sampled': true).
      -> $HADAMARD_WORK/cp_sweep.jsonl

The estimates only rank the structures (they are unreliable for small symmetric structures); cp_lowe.py searches."""
import argparse, json, math, sys, time
from math import comb, gcd
import numpy as np
from pathlib import Path as _P; sys.path[:0] = [str(_P(__file__).resolve().parent)]
import workdir as W
import cyclo_product as CP
import cp_lowe as CL


def is_prime(p):
    return p > 1 and all(p % d for d in range(2, int(p ** 0.5) + 1))


def cyclic_subgroups(m):
    """One generator per cyclic subgroup of Z_m^*."""
    seen, out = set(), []
    for k in range(1, m):
        if gcd(k, m) != 1: continue
        K, x = {1}, k
        while x != 1: K.add(x); x = x * k % m
        K = frozenset(K | {k % m})
        if K in seen: continue
        seen.add(K); out.append(k)
    return out


def structures(v, skew, minm, maxm):
    for q in range(3, v + 1):
        if v % q or not is_prime(q): continue
        m = v // q
        if (m < 2 and q != v) or m % q == 0: continue           # m = 1: v prime, the cyclotomy of GF(v) itself
        for e in range(2, q - 1, 2):
            if (q - 1) % e: continue
            F = (q - 1) // e
            if F < 2 or (skew and F % 2 == 0): continue
            if 4 * F - 4 > 4 * m * e: continue                    # a spike would exceed the budget alone
            for k in (cyclic_subgroups(m) if m > 1 else [1]):
                try: S = CP.Structure(v, q, e, skew=skew, kgen=k)
                except AssertionError: continue
                R = len(S.kreps)
                # members of a non-skew block without spikes: C(e, e/2) per representative (spikes give fewer)
                mem = comb(e, e // 2) ** R
                if mem < minm or mem > maxm * 30: continue
                v2 = (e & -e)
                for t in range(e):
                    if len(S.K) * t % e: continue
                    if t and len(S.K) == 1: continue
                    if skew and t % v2: continue
                    yield q, m, e, F, k, t, S, mem
    # prime-power fibres GF(p^k), k >= 2 (cp_field.FieldStructure; the group is Z_m x Z_p^k)
    import cp_field as FF
    for p, kk, e, F in FF.fibre_structures(v, skew):
        q = p ** kk; m = v // q
        if 4 * F - 4 > 4 * m * e: continue
        for k in (cyclic_subgroups(m) if m > 1 else [1]):
            try: S = FF.FieldStructure(v, p, kk, e, skew=skew, kgen=k)
            except AssertionError: continue
            R = len(S.kreps)
            mem = comb(e, e // 2) ** R
            if mem < minm or mem > maxm * 30: continue
            v2 = (e & -e)
            for t in range(e):
                if len(S.K) * t % e: continue
                if t and len(S.K) == 1: continue
                if skew and t % v2: continue
                yield q, m, e, F, k, t, S, mem


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=('skew', 'plain')); ap.add_argument('vs')
    ap.add_argument('--L', type=float, default=1e7); ap.add_argument('--ncls', type=int, default=6)
    ap.add_argument('--minm', type=float, default=1e5); ap.add_argument('--maxm', type=float, default=3e10)
    ap.add_argument('--samples', type=int, default=60000); ap.add_argument('--maxsec', type=float, default=1800)
    a = ap.parse_args()
    skew = a.mode == 'skew'
    out = W.WORK / 'cp_sweep.jsonl'
    done = set()
    if out.exists():
        for l in open(out):
            r = json.loads(l); done.add((r['v'], r['q'], r['e'], r['k'], r['twist'], r['mode'], r.get('config')))
    for v in [int(x) for x in a.vs.split(',')]:
        for q, m, e, F, k, t, S, mem in structures(v, skew, a.minm, a.maxm):
            base = dict(v=v, mode=a.mode, q=q, m=m, e=e, F=F, k=k, K=len(S.K), twist=t, reps=len(S.kreps), mem0=round(math.log10(mem), 1))
            if (v, q, e, k, t, a.mode, None) in done or (v, q, e, k, t, a.mode, 0) in done: continue
            t0 = time.time()
            sampled = False
            try:
                quads = CL.stage1_cached(S, not skew)
            except (RuntimeError, MemoryError):
                try:                                              # too large to list: sample a few classes instead
                    quads = CL.stage1_sample(S, a.ncls, np.random.default_rng(1), tries=300000); sampled = True
                except (RuntimeError, MemoryError) as ex:
                    rec = dict(base, error=str(ex)[:80]); print(json.dumps(rec), flush=True)
                    with open(out, 'a') as f: f.write(json.dumps(rec) + '\n')
                    continue
            base['sampled'] = sampled
            if not quads:
                rec = dict(base, classes=0, sec=round(time.time() - t0)); print(json.dumps(rec), flush=True)
                with open(out, 'a') as f: f.write(json.dumps(rec) + '\n')
                continue
            O = CL.Orbits(S, t)
            rng = np.random.default_rng(1)
            picks = sorted(set(np.linspace(0, len(quads) - 1, min(a.ncls, len(quads))).astype(int).tolist()))
            for c in picks:
                if time.time() - t0 > a.maxsec: break
                try:
                    P = CL.prepare(S, O, quads[c], a.samples, rng)
                    if S.skew and len(P['Pm0']) == 0: raise RuntimeError('no skew block member')
                    est, rank, off = CL.estimate_in_lists(P, lambda kk: a.L)
                    sizes = [round(math.log10(max(1, x)), 1) for x in P['sizes']]
                    big = sorted(P['sizes'])
                    # k-tree work per solution with A the smallest list and B, C, D of L members (pairs AB + CD)
                    Lb = [min(a.L, x) for x in big[1:]]
                    work = math.log10(big[0] * Lb[0] + Lb[1] * Lb[2]) - est
                    rec = dict(base, classes=len(quads), config=c, y=quads[c], rows=[int(sum(yy)) for yy in quads[c]],
                               logE=round(P['est0'], 2), logE_lists=round(est, 2), rank=rank, off=round(off, 2), sizes=sizes,
                               logwork=round(work, 2), sec=round(time.time() - t0))
                except (RuntimeError, MemoryError, AssertionError, ValueError) as ex:
                    rec = dict(base, classes=len(quads), config=c, error=str(ex)[:80])
                print(json.dumps(rec), flush=True)
                with open(out, 'a') as f: f.write(json.dumps(rec) + '\n')


if __name__ == '__main__':
    main()

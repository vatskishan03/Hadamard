"""Sweep the compression-chain search (chain.py) over the composite v of Cati-Pasechnik Table 5.

No Gaussian estimate is used.
- Scope: every subgroup H of Z_v^* with -1 not in H and MINO..MAXO orbits, and every row-sum pattern.
- Order: structures run in order of size.
- Per v: discovery mode stops at the first verified hit.
- Per structure: the result is 'found', 'ruled out' (exact) or 'nothing found (k combinations truncated)'.
- Output: one JSON line per structure and pattern, to $HADAMARD_WORK/chain_sweep.jsonl.

  python3 chain_sweep.py MINO MAXO v ...
"""
import json, math, sys, time
import numpy as np
from pathlib import Path as _P; sys.path[:0] = [str(_P(__file__).resolve().parent)]
import workdir as W
import spence_lib as L
from gs_sweep_subgroups import subgroups
import chain


def primes(v):
    return [p for p in range(2, v + 1) if v % p == 0 and all(p % q for q in range(2, int(p ** 0.5) + 1))]


def patterns(st, v):
    lim = int(math.isqrt(4 * v)); out = []
    for b in range(1, lim + 1, 2):
        for c in range(b, lim + 1, 2):
            d2 = 4 * v - 1 - b * b - c * c
            if d2 < c * c: continue
            d = math.isqrt(d2)
            if d * d == d2 and all(st.pool(r) > 0 for r in (b, c, d)): out.append((b, c, d))
    return out


def levels(v, H):
    """d1: the prime divisor on which H has the fewest orbits (strongest first filter); d2: the cofactor's largest prime."""
    ps = [p for p in primes(v) if p < v]
    if not ps: return None, None
    def orbits(d): return len({min(j * h % d for h in H) for j in range(d)})
    ps.sort(key=lambda p: orbits(p) / p)
    d1 = ps[0]
    rest = [p for p in ps if p != d1]
    d2 = max(rest) if rest else (d1 if v % (d1 * d1) == 0 and d1 * d1 < v else 0)
    if d2 == d1: d2 = d1 * d1 if v % (d1 ** 3) == 0 else 0
    return d1, d2


if __name__ == '__main__':
    mino, maxo = int(sys.argv[1]), int(sys.argv[2])
    out = W.WORK / 'chain_sweep.jsonl'
    for v in map(int, sys.argv[3:]):
        if len(primes(v)) == 1 and v == primes(v)[0]: continue          # prime: no compression available
        jobs = []
        for Hs, gens in subgroups(v).items():
            if (v - 1) in Hs: continue
            O = len({min(t * h % v for h in Hs) for t in range(v)})
            if not (mino <= O <= maxo): continue
            G = L.AbGroup([v]); st = L.Structure(G, [G.perm_from_matrix(np.array([[g]])) for g in gens])
            rep = np.array([np.flatnonzero(st.orb == o)[0] for o in range(st.O)])
            if np.any(st.orb[G.neg()[rep]][1:] == np.arange(1, st.O)): continue
            for pat in patterns(st, v):
                jobs.append((O, list(gens), sorted(Hs), pat))
        jobs.sort(key=lambda j: j[0])
        hit = False
        for O, gens, Hs, pat in jobs:
            d1, d2 = levels(v, Hs)
            if d1 is None: continue
            t = time.time()
            try:
                r = chain.run(v, gens, pat, d1, d2, exhaustive=False, verbose=False)
            except Exception as e:  # record and continue
                r = dict(status=f'error: {e!r}')
            rec = dict(v=v, gens=gens, H=len(Hs), O=O, pattern=list(pat), d1=d1, d2=d2, status=r['status'],
                       branches=r.get('branches'), truncated=r.get('truncated'), seconds=round(time.time() - t))
            with open(out, 'a') as f: f.write(json.dumps(rec) + '\n')
            print(json.dumps(rec), flush=True)
            if r.get('found'): hit = True; break
        print(f"== v={v}: {'HIT' if hit else 'no hit'} ({len(jobs)} structure-patterns)", flush=True)

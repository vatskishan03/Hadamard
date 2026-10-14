"""Exact compression-chain sweep (chain.py, plain Goethals-Seidel mode) over the unresolved Hadamard orders 4n with
composite n.

For every subgroup H of Z_n^* with MINO..MAXO orbits and every row-sum pattern a^2 + b^2 + c^2 + d^2 = 4n, it
records 'found', 'ruled out' (exact) or 'nothing found (truncated)' to $HADAMARD_WORK/had_chain_sweep.jsonl.

  python3 had_chain_sweep.py MINO MAXO n ...
"""
import json, math, sys, time
import numpy as np
from pathlib import Path as _P; sys.path[:0] = [str(_P(__file__).resolve().parent)]
import workdir as W
import spence_lib as L
from gs_sweep_subgroups import subgroups
import chain
from chain_sweep import primes, levels


def patterns(st, n):
    lim = int(math.isqrt(4 * n)); out = []
    for a in range(1, lim + 1, 2):
        for b in range(a, lim + 1, 2):
            for c in range(b, lim + 1, 2):
                d2 = 4 * n - a * a - b * b - c * c
                if d2 < c * c: continue
                d = math.isqrt(d2)
                if d * d == d2 and all(st.pool(r) > 0 for r in (a, b, c, d)): out.append((a, b, c, d))
    return out


if __name__ == '__main__':
    mino, maxo = int(sys.argv[1]), int(sys.argv[2])
    out = W.WORK / 'had_chain_sweep.jsonl'
    for n in map(int, sys.argv[3:]):
        if n in primes(n): continue
        jobs = []
        for Hs, gens in subgroups(n).items():
            O = len({min(t * h % n for h in Hs) for t in range(n)})
            if not (mino <= O <= maxo): continue
            G = L.AbGroup([n]); st = L.Structure(G, [G.perm_from_matrix(np.array([[g]])) for g in gens])
            for pat in patterns(st, n): jobs.append((O, list(gens), sorted(Hs), pat))
        jobs.sort(key=lambda j: j[0])
        hit = False
        for O, gens, Hs, pat in jobs:
            d1, d2 = levels(n, Hs)
            t = time.time()
            try:
                r = chain.run(n, gens, pat, exhaustive=False, verbose=False, budget=2400)
            except Exception as e:
                r = dict(status=f'error: {e!r}')
            rec = dict(n=n, gens=gens, H=len(Hs), O=O, pattern=list(pat), d1=d1, d2=d2, status=r['status'],
                       branches=r.get('branches'), truncated=r.get('truncated'), seconds=round(time.time() - t))
            with open(out, 'a') as f: f.write(json.dumps(rec) + '\n')
            print(json.dumps(rec), flush=True)
            if r.get('found'): hit = True; break
        print(f"== n={n}: {'HIT' if hit else 'no hit'} ({len(jobs)} structure-patterns)", flush=True)

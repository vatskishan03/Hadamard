"""enumerate origin-defect shapes: K=Z_m (m odd), q prime = 3 mod 4 (gcd 1), H_K odd order, nu<=NUMAX, with an origin cost combination sum 4m."""
import sys, os, json, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from shapes_lib import *
from math import gcd
def combos(m, q):
    A0s = [(q + 1) // 2, (q - 3) // 2]
    Os = sorted({0, (q + 1) // 2, (q - 3) // 2, q + 1, q - 3, q - 1})
    out = []
    for a in A0s:
        for b in Os:
            for c in Os:
                if c < b: continue
                for d in Os:
                    if d < c: continue
                    if a + b + c + d == 4 * m: out.append((a, b, c, d))
    return out
def shapes(mmin, mmax, NUMAX, NUMIN=3):
    res = []
    for m in range(mmin, mmax + 1, 2):
        subs = None
        for q in primes34(3, 8 * m + 3):
            if gcd(m, q) != 1: continue
            cb = combos(m, q)
            if not cb: continue
            if subs is None: subs = subgroups_odd(m)
            for H, gens in subs.items():
                nu = len(orbits(m, H))
                if nu > NUMAX or nu < NUMIN: continue
                res.append(dict(m=m, q=q, gens=gens, hsize=len(H), nu=nu, v=m * q, combos=cb))
    return res
if __name__ == '__main__':
    mmin, mmax, NUMAX = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    r = shapes(mmin, mmax, NUMAX)
    r.sort(key=lambda s: (s['nu'], s['v']))
    for s in r:
        print('nu', s['nu'], 'm', s['m'], 'q', s['q'], 'v', s['v'], 'H', s['gens'], '|H|', s['hsize'], 'ncombos', len(s['combos']), s['combos'][:2])
    print(len(r), file=sys.stderr)

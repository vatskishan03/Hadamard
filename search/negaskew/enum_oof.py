"""Origin-defect family of the quaternary Paley plug-in (qp.c): shapes (K = Z_m, q, H_K) of the Table 5 targets whose defect
costs can all sit at k = 0, i.e. 4|K| = A0 + B0 + C0 + D0 with A0 in {(q+1)/2, (q-3)/2} (the forced half defect of A) and B0, C0,
D0 in {0, (q+1)/2, (q-3)/2, q+1, q-3, q-1} (the pattern of both hits: every defect at the origin of K, NOTES.md section 89).
usage: python3 enum_oof.py          environment NUMAX (default 25: at most 25 K-orbits); PLUG_SETTLED=repo as in enum_plug.py
Writes runs/oof_shapes.json and prints one line per shape and the number of distinct v.  The targets that remain unsettled are
299, 413, 517 and 589; PLUG_SETTLED=all (all Table 5 entries m >= 4, v <= 999) also lists 209, 253, 341, 345, 391, 473 and 903, which
the repository settles (345 and 391 need NUMAX >= 15 and 17)."""
import json, itertools, sys, os
from math import gcd
from sympy import isprime
from enum_plug import T5, targets, subgroups_odd, orbits
NUMAX=int(os.environ.get('NUMAX','25'))
tg = sorted(v for v, m in T5.items() if m >= 4 and v <= 999) if os.environ.get('PLUG_SETTLED') == 'all' else targets()
res=[]
for v in tg:
    for q in range(3, v):
        if v % q or not isprime(q) or q % 4 != 3: continue
        m = v // q
        if m % 2 == 0 or m < 3 or gcd(m, q) != 1: continue
        A0s = [(q + 1) // 2, (q - 3) // 2]
        Os = sorted({0, (q + 1) // 2, (q - 3) // 2, q + 1, q - 3, q - 1})
        combos = []
        for a in A0s:
            for b in Os:
                for c in Os:
                    if c < b: continue
                    for d in Os:
                        if d < c: continue
                        if a + b + c + d == 4 * m: combos.append((a, b, c, d))
        if not combos: continue
        subs = subgroups_odd(m)
        for H, gens in sorted(subs.items(), key=lambda kv: -len(kv[0])):
            orbs = orbits(m, H); nu = len(orbs)
            if nu > NUMAX or nu < 3: continue
            res.append(dict(v=v, m5=T5[v], K=m, q=q, gens=gens, Hsize=len(H), nu=nu, combos=combos))
os.makedirs('runs', exist_ok=True)
json.dump(res, open('runs/oof_shapes.json', 'w'))
for s in sorted(res, key=lambda s: (s['nu'], s['v'])):
    print('nu', s['nu'], 'v', s['v'], 'm5', s['m5'], 'K', s['K'], 'q', s['q'], 'gens', s['gens'], '|H|', s['Hsize'], 'combos', s['combos'][:3], len(s['combos']))
print(len(res), 'shapes; distinct v:', len({s['v'] for s in res}))

"""Orbit-weight feasibility of the defect identity (EXACT, independent dynamic program).  The defect weights (e_p, e_m) of the identity (q+1) e_p + (q-3) e_m = 4|K| must be
realised by orbits of H_K: block A uses pairs {O, -O} (weight 2|O|), blocks B, C, D single orbits (weight |O|), every orbit at most once per block, type p or m.  If no target is
realisable the shape is EXACTLY empty ('identity-infeasible'); bbpat.enum_patterns must then also test 0 patterns (the program prints 'agrees' when the two verdicts coincide).
usage (from search/skewplug):  python3 patterns/check_weights.py                      the ten shapes of the non-cyclic Singer / Borel groups (1183, 1859, 2375, 1083, 1029, 2401)
                               python3 patterns/check_weights.py 'K' q 'hunits'       one shape, e.g.  "[37]" 23 "[[10]]"  (851: infeasible)
Expected for the default list: infeasible for '1859 Singer21' and '2375 S31' (the two EXACTLY empty cases of NOTES §88), feasible for the others."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
import json, warnings; warnings.filterwarnings('ignore')
from sk import SShape
from bbpat import identity_targets, enum_patterns
def reach(ws, maxe):
    S = {(0, 0)}
    for x in ws:
        T = set(S)
        for (a, b) in S:
            if a + x <= maxe: T.add((a + x, b))
            if b + x <= maxe: T.add((a, b + x))
        S = T
    return S
def feasible(sh):
    tg = identity_targets(sh)
    w = [len(p) for p in sh.spos]; A_w = [2 * w[o] for (o, po) in sh.apairs]
    q = sh.q
    if q == 3:
        maxp = max(ep for ep, em in tg); maxm = 4 * sh.N
        def reach3(ws):
            S = {(0, 0)}
            for x in ws:
                T = set(S)
                for (a, b) in S:
                    if a + x <= maxp: T.add((a + x, b))
                    if b + x <= maxm: T.add((a, b + x))
                S = T
            return S
        comb = lambda R1, R2: {(a + c, min(b + d, maxm)) for (a, b) in R1 for (c, d) in R2 if a + c <= maxp}
        R = comb(comb(comb(reach3(A_w), reach3(w)), reach3(w)), reach3(w))
        return any((ep, b) in R for (ep, em) in tg for b in range(maxm + 1)), len(tg)
    maxe = max(ep + em for ep, em in tg)
    RA = reach(A_w, maxe); RF = reach(w, maxe)
    comb = lambda R1, R2: {(a + c, b + d) for (a, b) in R1 for (c, d) in R2 if a + c <= maxe and b + d <= maxe}
    R = comb(comb(comb(RA, RF), RF), RF)
    return any(t in R for t in tg), len(tg)
S21 = [[11, 9], [2, 6]]; B39 = [[3, 3], [0, 3]]; S31 = [[1, 2, 3], [3, 1, 2], [4, 1, 3]]; S45 = [[12, 10], [14, 14]]
S57 = [[5, 0, 4], [6, 5, 0], [0, 6, 1]]; S171 = [[0, 3, 5], [0, 0, 3], [1, 4, 2]]
SHAPES = [([13, 13], 7, [S21], '1183 Singer21'), ([13, 13], 11, [S21], '1859 Singer21'), ([13, 13], 7, [B39], '1183 Borel39'), ([13, 13], 11, [B39], '1859 Borel39'),
                         ([5, 5, 5], 19, [S31], '2375 S31'), ([19, 19], 3, [S45], '1083 S45 (q=3)'), ([7, 7, 7], 3, [S57], '1029 S57'), ([7, 7, 7], 7, [S57], '2401 S57'), ([7, 7, 7], 3, [S171], '1029 S171'), ([7, 7, 7], 7, [S171], '2401 S171')]
if len(sys.argv) > 3: SHAPES = [(json.loads(sys.argv[1]), int(sys.argv[2]), json.loads(sys.argv[3]), sys.argv[1] + ' q=' + sys.argv[2] + ' ' + sys.argv[3])]
for (K, q, hu, name) in SHAPES:
    sh = SShape(K, q, hu)
    f, nt = feasible(sh)
    ntest = enum_patterns(sh, signs=False, maxleaves=2000, tlimit=60)['ntest']      # the pattern enumerator of bbpat.py: 0 patterns tested = identity-infeasible
    print(name, 'nu', sh.kCD.nu, 'weights', sorted(set(len(p) for p in sh.spos)), 'targets', nt, '-> weights realise a target:', f,
          '| bbpat tests', ntest, 'patterns:', 'agrees' if (ntest > 0) == f else 'DISAGREES', flush=True)

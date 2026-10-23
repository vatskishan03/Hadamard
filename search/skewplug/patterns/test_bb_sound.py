"""Soundness of the big-big filter (bbpat) on actual solutions: for every shape the exact allocation pipeline (../sk_exh.py) produces K-level solutions (up to 400 per shape);
each must (a) have an ALIVE position pattern under bbpat.check_pattern (with the count test whenever it is valid) and (b) satisfy the bound |sum_X BB_X(k)| <= N(k) at every lag
(direct_check).  Also provides pattern_of(ctx, sol), the pattern (orbit, type, sign) of a K-level solution (used by test_rowfilter.py and posctl_1100.py).
usage (from search/skewplug):  python3 patterns/test_bb_sound.py [i ...]       shapes by index (default all 13; the output line of each shape ends 'pattern NOT alive: 0 direct bound holds: n of n')"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
import json, warnings
warnings.filterwarnings('ignore')
import numpy as np
from sk import SShape
from sk_exh import exhaust
from bbpat import PatCtx, check_pattern

def pattern_of(ctx, sol):
    sh = ctx.sh; q = ctx.q; blocks = []
    for X, (s, t, e) in enumerate(sol):
        s = np.asarray(s); e = np.asarray(e); items = []
        if X == 0:
            for (o, po) in sh.apairs:
                p0 = int(ctx.pos[o][0])
                if e[p0] != 0: items.append((o, 0 if s[p0] == 0 else 1, +1)); items.append((po, 0 if s[p0] == 0 else 1, -1))
        else:
            for o in range(ctx.nu):
                p0 = int(ctx.pos[o][0])
                if e[p0] != 0: items.append((o, 0 if s[p0] == 0 else 1, +1))
        blocks.append(items)
    return blocks

def direct_check(sh, sol):
    N, q = sh.N, sh.q; bb = np.zeros(N, dtype=np.int64); Nk = np.zeros(N, dtype=np.int64)
    for (s, t, e) in sol:
        s = np.asarray(s, dtype=np.int64); e = np.asarray(e, dtype=np.int64); u = s + q * e
        D = (e != 0)
        for k in range(1, N):
            sk = sh.shift[k]            # j + k
            both = D & D[sk]
            bb[k] += (u * u[sk] * both).sum()
            # bounds
            Nk[k] += (np.abs(u) * (D & ~D[sk])).sum() + (np.abs(u[sk]) * (D[sk] & ~D)).sum() + ((~D) & (~D[sk])).sum()
    return bool((np.abs(bb[1:]) <= Nk[1:]).all())


SHAPES = [([3], 7, []), ([3], 11, []), ([5], 7, []), ([13], 7, [[3]]), ([5], 11, []), ([3], 19, []), ([5], 19, []), ([5], 23, []), ([5], 31, []), ([3], 23, []),
          ([7], 23, [[2]]), ([13], 19, [[3]]), ([13], 11, [[3]])]


if __name__ == '__main__':
    sel = [int(x) for x in sys.argv[1:]] or range(len(SHAPES))
    for i in sel:
        K, q, hu = SHAPES[i]
        sh = SShape(K, q, hu); ctx = PatCtx(sh)
        tot, sols = exhaust(sh, tmax=120, want_sols=400)
        nbad = 0; ndirect = 0; nrig = 0
        for sol in sols:
            blocks = pattern_of(ctx, sol)
            Et = sum(len(ctx.pos[o]) for b in blocks for (o, ty, sg) in b)
            rig = (q - 2) ** 2 > 2 * Et * q + 4 * sh.N
            nrig += int(rig)
            ok1, ok2, w = check_pattern(ctx, blocks, rigid=rig)
            if not (ok1 and ok2): nbad += 1
            ndirect += int(direct_check(sh, sol))
        print(i, K, q, hu, 'E_G', tot['E_G'], 'solutions checked', len(sols), 'pattern NOT alive:', nbad, 'direct bound holds:', ndirect, 'rigid-count-test applied:', nrig, flush=True)

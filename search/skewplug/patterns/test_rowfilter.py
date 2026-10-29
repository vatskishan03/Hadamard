"""Soundness of the exact row-sum pattern filter (sk_pat.row_patterns_skew): for actual K-level solutions (up to 300 per shape from the exact pipeline) the (allocation, m-split) class of the
solution must have at least one exact row-sum pattern.  Prints, for nine shapes, the number of solutions and the number of classes rejected (always 0: 1,900 solutions in the recorded run).
usage (from search/skewplug):  python3 patterns/test_rowfilter.py [i ...]       shapes by index 0..8"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
import warnings; warnings.filterwarnings('ignore')
import numpy as np
from sk import SShape
from sk_exh import exhaust
from bbpat import PatCtx, weights_ms
from sk_pat import row_patterns_skew
from test_bb_sound import pattern_of
SHAPES = [([3], 7, []), ([3], 11, []), ([5], 7, []), ([13], 7, [[3]]), ([5], 11, []), ([5], 19, []), ([5], 23, []), ([3, 3, 3], 7, [[[0, 0, 1], [0, 1, 1], [2, 1, 1]]]), ([7], 27, [[2]])]
sel = [int(x) for x in sys.argv[1:]] or range(len(SHAPES))
for (K, q, hu) in [SHAPES[i] for i in sel]:
    sh = SShape(K, q, hu); ctx = PatCtx(sh)
    tot, sols = exhaust(sh, tmax=100, want_sols=300)
    bad = 0
    for sol in sols:
        blocks = pattern_of(ctx, sol)
        A = tuple(sorted({(next(j for j, (o, po) in enumerate(sh.apairs) if it[0] in (o, po)), it[1]) for it in blocks[0]}))
        free = [tuple(sorted((it[0], it[1]) for it in b)) for b in blocks[1:]]
        al, ms = weights_ms(sh, A, free[0], free[1], free[2])
        if not row_patterns_skew(sh, al, ms): bad += 1
    print(K, q, hu, 'solutions', len(sols), 'classes rejected by the row filter:', bad, flush=True)

"""POSITIVE CONTROL of the pattern pipeline for matrix-valued H_K and non-cyclic K: the repository's skew-Hadamard matrix skew/1100 (K = Z5 x Z5, q = 11, H_K = <Singer matrix of order 3>,
S = squares of F_11; v = 275) is converted to K-level blocks (s, t, e) and checked exactly at K level and at G level; its defect pattern must be alive under the big-big filter, and the
exhaustive per-pattern pipeline (explicit lists + exact join + exact level 2) must RECOVER it.  Expected output: both checks True, 'bound ... holds: True', 'signed-pattern alive: (True, True)',
the pipeline record with n1 = 64 (level-1 instances) and compl = 40 (completions), and 'recovered solution K-level check: True' for each recovered solution (about 10 s).
usage (from search/skewplug):  python3 patterns/posctl_1100.py"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
import json, warnings, numpy as np
warnings.filterwarnings('ignore')
from sk import SShape
from pat_exh import solve_pattern
from gl1 import L1
from sl2io import header_text
from bbpat import PatCtx, check_pattern
c = json.load(open(os.path.join(REPO, 'skew', '1100', 'certificate.json')))
seqs = [np.array([1 if ch == '+' else -1 for ch in s]) for s in c['sequences']]
sh = SShape([5, 5], 11, [[[2, 4], [2, 2]]])
N, q = sh.N, 11
blocks = []
for s_ in seqs:
    x0 = s_[0::q]; x1 = s_[1::q]; x2 = s_[2::q]
    # H-invariance check: x1 constant on squares, x2 on non-squares
    sq = [1, 3, 4, 5, 9]; nsq = [2, 6, 7, 8, 10]
    for k in range(N):
        assert all(s_[k * q + f] == x1[k] for f in sq) and all(s_[k * q + f] == x2[k] for f in nsq)
    e = (x1 + x2) // 2; t = (x1 - x2) // 2; s = x0 - e
    blocks.append((s, t, e))
print('K-level check (exact, all k):', sh.check_solution(blocks), ' G-level skew+target check:', sh.g_target_skew_ok(sh.glevel_blocks(blocks)))
ctx = PatCtx(sh)
from test_bb_sound import pattern_of, direct_check
pat = pattern_of(ctx, blocks)
print('pattern (orbit, type, sign):', pat, ' weights per block:', [sum(len(ctx.pos[o]) for (o, ty, sg) in b) for b in pat])
print('bound |sum BB| <= N(k) holds:', direct_check(sh, blocks), ' signed-pattern alive:', check_pattern(ctx, pat)[:2])
# pattern -> pipeline
A = tuple(sorted({(next(j for j, (o, po) in enumerate(sh.apairs) if it[0] in (o, po)), it[1]) for it in pat[0]}))
free = [tuple(sorted((it[0], it[1]) for it in b)) for b in pat[1:]]
l1 = L1(sh); hdr = header_text(sh)
rec, sols = solve_pattern(sh, l1, hdr, A, free[0], free[1], free[2], want_sols=3)
print('pipeline on the pattern:', rec)
for s in sols: print('recovered solution K-level check:', sh.check_solution(s))

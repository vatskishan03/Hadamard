"""Soundness of the big-big filter on a rigid shape with NONZERO level 1: the sum over all alive patterns (no symmetry reduction, no row-sum filter) of the ordered level-1 instances, times 8,
equals the allocation pipeline's N1 (exact).  usage (from search/skewplug):  python3 patterns/val_n1.py 'K' q 'hunits' EXPECTED_N1
Expected: "[33]" 23 "[[4]]" 124416 (759: Z33 with <4> x F23; the filter removes 180 of 486 patterns and no level-1 instance)."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
import json, warnings, math
warnings.filterwarnings('ignore')
from collections import Counter
from sk import SShape
from bbpat import enum_patterns
from pat_exh import solve_pattern
from gl1 import L1
from sl2io import header_text
K = json.loads(sys.argv[1]); q = int(sys.argv[2]); hu = json.loads(sys.argv[3]); expect = int(sys.argv[4])
sh = SShape(K, q, hu); l1 = L1(sh); hdr = header_text(sh)
r = enum_patterns(sh)
print('alive patterns', len(r['alive']), 'killed(signed)', r['kill_sign'], 'killed(count)', r['kill_count'], flush=True)
tot = 0
for a in r['alive']:
    rec, sols = solve_pattern(sh, l1, hdr, a['A'], a['B'], a['C'], a['D'], want_sols=0, tl2=1.0)
    mult = 6; c = Counter([a['B'], a['C'], a['D']])
    for v in c.values(): mult //= math.factorial(v)
    tot += rec['n1'] * mult
print('ordered level-1 instances x 8 =', 8 * tot, 'expected N1_G', expect, 'MATCH' if 8 * tot == expect else 'MISMATCH', flush=True)

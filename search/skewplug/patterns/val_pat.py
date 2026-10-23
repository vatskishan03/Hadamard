"""Validation of the pattern pipeline: the sum over ALL alive position patterns (no symmetry reduction, no row-sum filter) of the ordered level-2 completions, times 8, equals the exact count
of the allocation pipeline (../sk_exh.py).  usage (from search/skewplug):  python3 patterns/val_pat.py 'K' q 'hunits' EXPECTED
Expected counts: "[13]" 7 "[[3]]" 930816 (Z13 with <3> x F7);  "[3,3,3]" 7 "[[[0,0,1],[0,1,1],[2,1,1]]]" 6144 (Z3^3 with a Singer matrix of order 13 x F7);  "[7]" 27 "[[2]]" 6144 (Z7 x GF(27)).
Last line: 'ordered completions x 8 = N expected N MATCH'."""
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
r = enum_patterns(sh, signs=False)
print('patterns', len(r['alive']), flush=True)
tot = 0; n1 = 0
for a in r['alive']:
    rec, sols = solve_pattern(sh, l1, hdr, a['A'], a['B'], a['C'], a['D'], want_sols=0)
    mult = math.factorial(3)
    c = Counter([a['B'], a['C'], a['D']])
    for v in c.values(): mult //= math.factorial(v)
    tot += rec['compl'] * mult; n1 += rec['n1'] * mult
print('ordered completions x 8 =', 8 * tot, 'expected', expect, 'MATCH' if 8 * tot == expect else 'MISMATCH', flush=True)

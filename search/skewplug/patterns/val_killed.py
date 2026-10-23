"""Empirical soundness of the big-big filter on a shape it declares EXACTLY EMPTY: the exhaustive per-pattern pipeline (pat_exh) is run on a random sample of the KILLED patterns (the count test
of the rigid regime is kept, the signed test is switched off); every one must have 0 level-1 instances.
usage (from search/skewplug):  python3 patterns/val_killed.py 'K' q 'hunits' NSAMPLE      e.g.  "[39]" 23 "[[16]]" 120  (897: 120 of the 7,332 killed patterns, 0 level-1 instances in all)"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
import json, warnings, random
warnings.filterwarnings('ignore')
from sk import SShape
from bbpat import enum_patterns
from pat_exh import solve_pattern
from gl1 import L1
from sl2io import header_text
K = json.loads(sys.argv[1]); q = int(sys.argv[2]); hu = json.loads(sys.argv[3]); nsamp = int(sys.argv[4])
sh = SShape(K, q, hu); l1 = L1(sh); hdr = header_text(sh)
r = enum_patterns(sh, signs=False)
pats = r['alive']; random.Random(1).shuffle(pats)
print('patterns without the signed test:', len(pats), 'count-killed:', r['kill_count'], flush=True)
bad = 0; tot_n1 = 0
for a in pats[:nsamp]:
    rec, _ = solve_pattern(sh, l1, hdr, a['A'], a['B'], a['C'], a['D'], want_sols=0, tl2=1.0)
    tot_n1 += rec['n1']; bad += int(rec['n1'] > 0)
print(K, q, hu, 'sample', min(nsamp, len(pats)), 'killed patterns with level-1 instances (must be 0):', bad, 'total instances', tot_n1, flush=True)

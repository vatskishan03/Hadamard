"""Validation of the symmetry reduction of patterns (patsym.Sym): over ALL alive patterns of a shape the pattern and its canonical representative have the same zero-ness of the exhaustive
per-pattern result (completions > 0).  Prints the number of patterns, of canonical classes, of classes with nonzero completions and of INCONSISTENT classes (always 0).
usage (from search/skewplug):  python3 patterns/val_sym.py 'K' q 'hunits'      e.g.  "[3]" 7 "[]";  "[3]" 11 "[]";  "[5]" 7 "[]";  "[13]" 7 "[[3]]" """
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
import json, warnings
warnings.filterwarnings('ignore')
from collections import defaultdict
from sk import SShape
from bbpat import enum_patterns
from patsym import Sym
from pat_exh import solve_pattern
from gl1 import L1
from sl2io import header_text
K = json.loads(sys.argv[1]); q = int(sys.argv[2]); hu = json.loads(sys.argv[3])
sh = SShape(K, q, hu); l1 = L1(sh); hdr = header_text(sh); sym = Sym(sh)
r = enum_patterns(sh, signs=False)
groups = defaultdict(list)
for a in r['alive']:
    key = sym.canon(a['A'], a['B'], a['C'], a['D'])
    rec, _ = solve_pattern(sh, l1, hdr, a['A'], a['B'], a['C'], a['D'], want_sols=0)
    groups[key].append(rec['compl'] > 0)
bad = sum(1 for k, v in groups.items() if len(set(v)) > 1)
nz = sum(1 for k, v in groups.items() if any(v))
print(K, q, hu, 'patterns', len(r['alive']), 'canonical classes', len(groups), 'classes with nonzero completions', nz, 'classes with INCONSISTENT zero-ness', bad, 'fixed pts', len(sym.fp), 'units', len(sym.umaps), flush=True)

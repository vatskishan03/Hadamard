"""Big-big filter + exact row-sum filter + symmetry reduction of the defect-position patterns of a shape, without any join (no level 1).
Prints one RESULT line; 'after row filter 0' means that every pattern dies: the shape is EXACTLY empty ('tested 0' means that the orbit weights cannot realise
the defect identity at all).  usage (from search/skewplug):  python3 patterns/run_bb.py 'K' q 'hunits' [tlimit]     e.g.  python3 patterns/run_bb.py "[39]" 23 "[[16]]" 900"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
import json, time, warnings
warnings.filterwarnings('ignore')
from sk import SShape
from bbpat import enum_patterns, rowfilter
from patsym import Sym
K = json.loads(sys.argv[1]); q = int(sys.argv[2]); hu = json.loads(sys.argv[3])
sh = SShape(K, q, hu)
r = enum_patterns(sh, tlimit=float(sys.argv[4]) if len(sys.argv) > 4 else 900)
alive, removed = rowfilter(sh, r['alive'])
r['alive'] = alive
sym = Sym(sh); red = sym.reduce(alive)
print('RESULT', K, q, hu, 'nu', sh.kCD.nu, 'targets', r['targets'], 'tested', r['ntest'], 'killed(count)', r['kill_count'], 'killed(signed)', r['kill_sign'], 'alive', len(alive) + removed, 'after row filter', len(alive), 'mod symmetry', len(red), 'fixed pts', len(sym.fp), 'units', len(sym.umaps), 'timed_out', r['timed_out'], 'secs', round(r['secs'], 1), flush=True)

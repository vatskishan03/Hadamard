"""tokens_est4.py -- est4 (HEURISTIC) log10 E and log10 N1 of every tau-feasible allocation token of a Spence plug-in shape, and its t bits (README, "Order 5912").
usage: tokens_est4.py Kmods r hg OUT.json      e.g. tokens_est4.py '[211]' 7 '[[[19],1]]' ../data/gs/p5912/tokenE_739.json"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gs_samp import make_shape
from prod_tau import tokens_for
from est4 import Est4
from k4 import Shape4

Kmods = json.loads(sys.argv[1]); r = int(sys.argv[2]); hg = json.loads(sys.argv[3]); out = sys.argv[4]
t0 = time.time(); toks = tokens_for(make_shape('spence', Kmods, r, hg, None))
print('tokens', len(toks), '%.1fs' % (time.time() - t0), flush=True)
es = Est4(Shape4(Kmods, r, [(tuple(u), ph) for u, ph in hg], ([], [()])), M=1000, seed=1)
res = []
for al in toks:
    a = tuple(sorted(n for n, m in al))
    try: le, _, tb = es.single('all', a); n1, _, _ = es.single('u', a)
    except Exception: le, n1, tb = None, None, None
    res.append(dict(token=[list(x) for x in al], alloc=list(a), logE=None if le is None else float(le), logN1=None if n1 is None else float(n1), tbits=None if tb is None else float(tb)))
    print(a, 'logE', le, 'logN1', n1, 'tbits', tb, flush=True)
json.dump(res, open(out, 'w'))

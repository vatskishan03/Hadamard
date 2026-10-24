"""census_a.py -- enumerate Paley plug-in shapes (modes a/b/c of the real plug-in shape census of search/NOTES.md section 86, with rmin = 3 available) for n = 2p and apply the exact identity / row-sum filters.
usage: census_a.py N OUT.jsonl [Omin Omax] [product [qkinds [modes [rmin]]]]   (the commands of the census of search/NOTES.md section 92 are in README.md)"""
import sys, json, time
import os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'gs'))
from shapes4 import *
import shapes4b

def run(n, out, Omin=2, Omax=28, product=False, qkinds=('z3','z2z2'), modes='abc', rmin=5):
    f = open(out, 'a'); t0 = time.time(); cnt = 0; live = 0
    gen = shapes4b.enumerate_shapes_b(n, modes=modes, Orange=(Omin, Omax), qkinds=qkinds, rmin=rmin) if product else enumerate_shapes(n, modes=modes, Orange=(Omin, Omax), qkinds=qkinds, rmin=rmin)
    for d in gen:
        hg = [(tuple(u), ph) for u, ph in d['hg']]
        try:
            sh = Shape4(d['K'], d['r'], hg, d['qspec'])
        except AssertionError:
            d['status'] = 'killed:phi-not-invariant'; f.write(json.dumps(d) + '\n'); continue
        d['Did'] = int(sh.Did); d['rowsq'] = int(sh.rowsq); d['nU'] = sh.nU; d['R'] = sh.R
        lo, hi = (-(-sh.Did // (sh.q + 1)), (sh.Did // (sh.q - 3) if sh.q > 3 else 4 * sh.N))
        d['defwin'] = [lo, hi]
        if sh.Did < 0: d['status'] = 'killed:identity-negative'
        elif sh.kCD.fw > hi: d['status'] = 'killed:forced-defects-exceed'
        else:
            ok, why = dp_filter(sh)
            d['status'] = 'live' if ok else 'killed:' + why
        live += d['status'] == 'live'; cnt += 1
        f.write(json.dumps(d, default=lambda o: o.item() if hasattr(o, 'item') else str(o)) + '\n'); f.flush()
    print(n, 'product' if product else 'cyclic', 'shapes', cnt, 'live', live, '%.0fs' % (time.time() - t0), flush=True)

if __name__ == '__main__':
    n = int(sys.argv[1]); out = sys.argv[2]
    Omin = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    Omax = int(sys.argv[4]) if len(sys.argv) > 4 else 28
    prod = len(sys.argv) > 5 and sys.argv[5] == 'product'
    qk = tuple(sys.argv[6].split(',')) if len(sys.argv) > 6 else ('z3', 'z2z2')
    modes = sys.argv[7] if len(sys.argv) > 7 else 'abc'
    rmin = int(sys.argv[8]) if len(sys.argv) > 8 else 5
    run(n, out, Omin, Omax, prod, qk, modes, rmin)

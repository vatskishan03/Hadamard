"""calib_odd.py -- calibration of the Gaussian/lattice estimator est4 on the exact counts of odd-|G| Spence shapes (val_odd.jsonl: K-level engine == G-level brute force)."""
import sys, json, math, signal
import numpy as np
import os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import k4, est4
class TO(Exception): pass
def _al(a, b): raise TO()
signal.signal(signal.SIGALRM, _al)
res = []
for l in open(os.path.join(HERE, '..', 'data', 'gs', 'val_odd.jsonl')):
    d = json.loads(l)
    if d['status'] != 'MATCH': continue
    hg = [(tuple(u), ph) for u, ph in d['hg']]
    try:
        sh = k4.Shape4(d['K'], d['r'], hg, ([], [()]))
        signal.alarm(60)
        try: e = est4.estimate_shape(sh, M=400, nfull=100, topk=4, ntop=1)
        finally: signal.alarm(0)
    except TO: continue
    except Exception as ex: continue
    le = e.get('logE')
    exact = d['G_count'] / 16.0
    res.append(dict(K=d['K'], r=d['r'], hg=d['hg'], est=le, exact=(math.log10(exact) if exact > 0 else None)))
nz = [r for r in res if r['exact'] is not None and r['est'] is not None and np.isfinite(r['est'])]
z = [r for r in res if r['exact'] is None and r['est'] is not None]
print('shapes', len(res), 'nonzero exact', len(nz), 'zero exact', len(z))
if nz:
    diff = np.array([r['est'] - r['exact'] for r in nz])
    print('nonzero: est - exact (decades): mean %.2f sd %.2f min %.2f max %.2f' % (diff.mean(), diff.std(), diff.min(), diff.max()))
fa = [r['est'] for r in z if np.isfinite(r['est']) and r['est'] > 0]
print('zero-count shapes with est logE > 0 (false alarms):', len(fa), 'of', len(z), ' max est', max(fa) if fa else None)
print('zero-count shapes with est logE finite:', sum(1 for r in z if np.isfinite(r['est'])))
json.dump(res, open('calib_odd.json', 'w'))

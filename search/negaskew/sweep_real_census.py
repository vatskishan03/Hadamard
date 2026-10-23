"""Real-alphabet sweep of the plug-in engine over the census-undecided shapes (census_shapes.py: odd v < 2500 without a matrix of order 4v,
K = Z_m, q prime = 3 mod 4, H_K odd order) with NUMIN <= nu <= NUMAX K-orbits: `qp real=1 l2c=1 first=1` per shape with a time limit (a skew
order 4v would be a census hit).  Resumable.
usage: python3 sweep_real_census.py     environment: NUMIN (7), NUMAX (10), TMO seconds per shape (40), LOGF (runs/real_census.jsonl)
  NUMIN=7 NUMAX=10 TMO=40 python3 sweep_real_census.py                       139 shapes with 7 and 9 orbits: 128 empty, 11 timeouts, no hit
  NUMIN=11 NUMAX=15 TMO=20 LOGF=runs/real_census_hi.jsonl python3 sweep_real_census.py     151 shapes: 65 empty, 86 timeouts, no hit
REV=1 runs the shapes in reverse order.  Prints HIT / TIMEOUT / empty per shape; a TIMEOUT is UNKNOWN and depends on the speed of the machine.
(The shapes with at most 5 orbits were swept by sweep.py with mode 0.)  Needs ./qp."""
import json, os, sys, subprocess, time
from qshape import kshape, write_shape
from census_shapes import shapes as census_shapes
NUMIN = int(os.environ.get('NUMIN', '7')); NUMAX = int(os.environ.get('NUMAX', '10')); TMO = int(os.environ.get('TMO', '40'))
LOG = os.environ.get('LOGF', 'runs/real_census.jsonl')
os.makedirs('runs', exist_ok=True); os.makedirs('shapes', exist_ok=True)
done = set()
if os.path.exists(LOG):
    for l in open(LOG):
        try: d = json.loads(l); done.add((d['v'], d['K'], d['q'], tuple(d['gens'])))
        except Exception: pass
shapes = [s for s in census_shapes(NUMAX) if s[0] >= NUMIN]
if os.environ.get('REV'): shapes = shapes[::-1]
print(len(shapes), 'shapes', file=sys.stderr)
for nu, v, m, q, gens, hs in shapes:
    key = (v, m, q, tuple(gens))
    if key in done: continue
    fn = f'shapes/c_{v}_{m}_{q}_{"_".join(map(str,gens))}.txt'
    sh = kshape([m], gens); write_shape(sh, q, fn)
    t0 = time.time()
    try:
        r = subprocess.run(['./qp', fn, 'real=1', 'l2c=1', 'first=1', 'maxprint=1'], capture_output=True, text=True, timeout=TMO)
        out, err, rc = r.stdout, r.stderr, r.returncode
    except subprocess.TimeoutExpired as e:
        out = e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or ''); err = e.stderr.decode() if isinstance(e.stderr, bytes) else (e.stderr or ''); rc = 'timeout'
    sol = [l for l in out.split('\n') if l.startswith('SOL')]
    stats = [l for l in out.split('\n') if l.startswith('STATS')]
    lists = [l for l in err.split('\n') if l.startswith('lists:')]
    rec = dict(v=v, K=m, q=q, gens=gens, Hsize=hs, nu=nu, rc=str(rc), sec=round(time.time() - t0, 1), sol=sol[0] if sol else None, stats=stats[0][:300] if stats else None, lists=lists[0] if lists else None)
    open(LOG, 'a').write(json.dumps(rec) + '\n')
    flag = 'HIT' if sol else ('TIMEOUT' if rc == 'timeout' else 'empty')
    print(flag, nu, v, m, q, gens, rec['sec'], (stats[0][:140] if stats else ''), flush=True)

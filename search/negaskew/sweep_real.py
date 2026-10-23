"""Real-alphabet sweep of the plug-in engine: for every Table 5 target v (m >= 4, v <= 999, 4v and 8v not in skew/; enum_plug.targets())
and every shape (K = Z_m, q prime = 3 mod 4, gcd(m, q) = 1, H_K an odd-order subgroup including the trivial one) with NUMIN <= nu <= NUMAX
K-orbits, run `qp real=1 l2c=1 first=1` (the real skew class of search/skewplug, skew order 4v; exhaustive, stops at the first solution) with a
time limit.  Resumable: the shapes already in the log are skipped.
usage: python3 sweep_real.py        environment: NUMIN (3), NUMAX (10), TMO seconds per shape (120), LOGF (runs/real_sweep.jsonl)
  NUMAX=10 TMO=120 python3 sweep_real.py                       85 shapes with 3..10 orbits: 81 empty, 4 timeouts (623, 627, 759, 959), no hit
  NUMIN=11 NUMAX=15 TMO=30 LOGF=runs/real_hi.jsonl python3 sweep_real.py     21 shapes: 9 empty (517, 649, 737, 869, 913, 299, 611, 589, 851), no hit
Prints HIT / TIMEOUT / empty per shape; a HIT line is a skew order 4v (certify with cert_real.py).  Needs ./qp (README.md); the timeouts
depend on the speed of the machine (a TIMEOUT is UNKNOWN, not a result).  PLUG_SETTLED=repo: see enum_plug.py."""
import json, os, sys, subprocess, time
from math import gcd
from sympy import isprime
from enum_plug import T5, targets, subgroups_odd, orbits
from qshape import kshape, write_shape
NUMAX = int(os.environ.get('NUMAX', '10')); NUMIN = int(os.environ.get('NUMIN', '3')); TMO = int(os.environ.get('TMO', '120'))
LOG = os.environ.get('LOGF', 'runs/real_sweep.jsonl')
os.makedirs('runs', exist_ok=True); os.makedirs('shapes', exist_ok=True)
done = set()
if os.path.exists(LOG):
    for l in open(LOG):
        try: d = json.loads(l); done.add((d['v'], d['K'], d['q'], tuple(d['gens'])))
        except Exception: pass
shapes = []
for v in targets():
    for q in range(3, v):
        if v % q or not isprime(q) or q % 4 != 3: continue
        m = v // q
        if m % 2 == 0 or m < 3 or gcd(m, q) != 1: continue
        subs = subgroups_odd(m)
        for H, gens in subs.items():
            nu = len(orbits(m, H))
            if nu < NUMIN or nu > NUMAX: continue
            shapes.append((nu, v, m, q, list(gens), len(H)))
shapes.sort()
print(len(shapes), 'shapes', file=sys.stderr)
for nu, v, m, q, gens, hs in shapes:
    key = (v, m, q, tuple(gens))
    if key in done: continue
    fn = f'shapes/r_{v}_{m}_{q}_{"_".join(map(str,gens))}.txt'
    sh = kshape([m], gens); write_shape(sh, q, fn)
    t0 = time.time()
    try:
        r = subprocess.run(['./qp', fn, 'real=1', 'l2c=1', 'first=1', 'maxprint=1'], capture_output=True, text=True, timeout=TMO)
        out, err, rc = r.stdout, r.stderr, r.returncode
    except subprocess.TimeoutExpired as e:
        out, err, rc = (e.stdout or b'').decode() if isinstance(e.stdout, bytes) else (e.stdout or ''), (e.stderr or b'').decode() if isinstance(e.stderr, bytes) else (e.stderr or ''), 'timeout'
    sol = [l for l in out.split('\n') if l.startswith('SOL')]
    stats = [l for l in out.split('\n') if l.startswith('STATS')]
    lists = [l for l in err.split('\n') if l.startswith('lists:')]
    rec = dict(v=v, m5=T5[v], K=m, q=q, gens=gens, Hsize=hs, nu=nu, rc=str(rc), sec=round(time.time() - t0, 1), sol=sol[0] if sol else None, stats=stats[0][:300] if stats else None, lists=lists[0] if lists else None)
    open(LOG, 'a').write(json.dumps(rec) + '\n')
    flag = 'HIT' if sol else ('TIMEOUT' if rc == 'timeout' else 'empty')
    print(flag, nu, v, m, q, gens, rec['sec'], (stats[0][:140] if stats else ''), flush=True)

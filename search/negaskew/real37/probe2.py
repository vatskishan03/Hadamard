"""level-2 timing probe (binary: $QP, default ./qp37r): python3 probe2.py JOBS.json OUT.jsonl L2P TMO  -> per job: l1 instances, sampled level-2 instances, t_l2, estimated full level-2 time, l2_toobig"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)      # ../qshape.py is the base shape code
import json, os, sys, subprocess, time, resource, re
from qshape import kshape
from qshape37 import write_shape
jobs = json.load(open(sys.argv[1])); LOG = sys.argv[2]; P = float(sys.argv[3]); TMO = int(sys.argv[4]); os.makedirs('shapes', exist_ok=True)
done = set()
if os.path.exists(LOG):
    for l in open(LOG):
        try: done.add(json.loads(l)['name'])
        except Exception: pass
def lim():
    resource.setrlimit(resource.RLIMIT_AS, (3500000 * 1024, 3500000 * 1024)); os.nice(12)
for j in jobs:
    if j['name'] in done: continue
    fn = f'shapes/{j["name"]}.txt'
    if not os.path.exists(fn):
        sh = kshape(j.get('Kmods') or [j['K']], j['gens']); write_shape(sh, j['q'], fn)
    t0 = time.time()
    try:
        r = subprocess.run([os.environ.get('QP', './qp37r'), fn, 'real=1', 'l2c=1', 'sym=1', 'uA=1', 'rB=1', 's3=1', 'first=1', 'maxprint=1', f'l2p={P}', 'seed=5'], capture_output=True, text=True, timeout=TMO, preexec_fn=lim); out, err, rc = r.stdout, r.stderr, r.returncode
    except subprocess.TimeoutExpired as e:
        out = e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or ''); err = e.stderr.decode() if isinstance(e.stderr, bytes) else (e.stderr or ''); rc = 'timeout'
    st = [l for l in out.split('\n') if l.startswith('STATS')]; rec = dict(name=j['name'], rc=str(rc), sec=round(time.time() - t0, 1), stats=st[0] if st else None)
    if st:
        g = lambda k: float(re.search(k + r'=([\d.]+)', st[0]).group(1))
        rec.update(l1=g('l1_instances'), l2inst=g('l2_instances'), t_l2=g('t_l2'), toobig=g('l2_toobig'), t_total=g('time'), est_l2_full=round(g('t_l2') / P, 0), est_total=round(g('time') + g('t_l2') / P - g('t_l2'), 0), l2sol=g('l2_solutions'))
    else:
        pr = [l for l in err.split('\n') if l.startswith('PROG ')]; rec['lastprog'] = pr[-1] if pr else None
    open(LOG, 'a').write(json.dumps(rec) + '\n')
    print(j['name'], rec.get('rc'), rec['sec'], {k: rec.get(k) for k in ('l1', 'l2inst', 't_l2', 'toobig', 'est_total')}, flush=True)

"""dry runs (exact level-1 probe totals) with the symmetry reductions (uA tF sym s3) for a job list: python3 dry2.py JOBS.json OUT.jsonl TMO   (binary: $QP, default ./qp37r)"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)      # ../qshape.py is the base shape code
import json, os, sys, subprocess, time, resource
from qshape import kshape
from qshape37 import write_shape
jobs = json.load(open(sys.argv[1])); LOG = sys.argv[2]; TMO = int(sys.argv[3]); os.makedirs('shapes', exist_ok=True)
done = set()
if os.path.exists(LOG):
    for l in open(LOG):
        try: done.add(json.loads(l)['name'])
        except Exception: pass
def lim():
    resource.setrlimit(resource.RLIMIT_AS, (3500000 * 1024, 3500000 * 1024)); os.nice(8)
for j in jobs:
    if j['name'] in done: continue
    fn = f'shapes/{j["name"]}.txt'
    if not os.path.exists(fn):
        sh = kshape(j.get('Kmods') or [j['K']], j['gens']); write_shape(sh, j['q'], fn)
    t0 = time.time()
    try:
        r = subprocess.run([os.environ.get('QP', './qp37r'), fn, 'real=1', 'l2c=1', 'dry=1', 'sym=1', 'uA=1', 'tF=1', 's3=1'], capture_output=True, text=True, timeout=TMO, preexec_fn=lim); out, err, rc = r.stdout, r.stderr, r.returncode
    except subprocess.TimeoutExpired as e:
        out = e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or ''); err = e.stderr.decode() if isinstance(e.stderr, bytes) else (e.stderr or ''); rc = 'timeout'
    dry = [l for l in out.split('\n') if l.startswith('DRY')]; lists = [l for l in err.split('\n') if l.startswith('lists')]
    rec = dict(name=j['name'], v=j['v'], K=j['K'], q=j['q'], gens=j['gens'], nu=j.get('nu'), rc=str(rc), sec=round(time.time() - t0, 1), dry=dry[0] if dry else None, lists=lists[0] if lists else None, errtail=err.strip().split('\n')[-1][:160])
    open(LOG, 'a').write(json.dumps(rec) + '\n'); print(j['nu'], j['v'], j['K'], j['q'], j['gens'], rc, rec['sec'], dry[0] if dry else rec['errtail'], flush=True)

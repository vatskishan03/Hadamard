"""exhaustive real-mode queue (search/NOTES.md section 90).  usage: python3 runq.py JOBS.json LOG.jsonl
JOBS.json: list of dicts {name, v, K, q, gens, tmo, extra (list of 'key=val'), [Kmods]}.  Each job runs $QP (default ./qp37r) shape real=1 l2c=1 first=1 maxprint=1 sym=1 uA=1 tF=1 s3=1 prog=60 EXTRA under ulimit -v 3500000, nice 5, timeout tmo (s).
Resume-safe: a job whose name has a record with status in (EMPTY, HIT, INCOMPLETE, ERROR) in LOG is skipped; TIMEOUT records are retried only when the job's tmo is larger than the recorded one.
Status: HIT (a SOL line; hits/NAME.sol written), EMPTY (rc 0, no SOL, l2_toobig = 0: EXHAUSTIVE within the class, symmetry-reduced by C<->D), INCOMPLETE (rc 0 but l2_toobig > 0), TIMEOUT, ERROR."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)      # ../qshape.py is the base shape code
import json, os, sys, subprocess, time, signal, re
from qshape import kshape
from qshape37 import write_shape
jobs = json.load(open(sys.argv[1])); LOG = sys.argv[2]; QP = os.environ.get('QP', './qp37r')
for d_ in ('shapes', 'runs', 'logs', 'hits'): os.makedirs(d_, exist_ok=True)
def load():
    d = {}
    if os.path.exists(LOG):
        for l in open(LOG):
            try: r = json.loads(l); d[r['name']] = r
            except Exception: pass
    return d
for j in jobs:
    prev = load().get(j['name'])
    if prev and (prev['status'] in ('EMPTY', 'HIT', 'INCOMPLETE', 'ERROR') or (prev['status'] == 'TIMEOUT' and prev['tmo'] >= j['tmo'])): continue
    Km = j.get('Kmods') or [j['K']]
    fn = f'shapes/{j["name"]}.txt'
    if not os.path.exists(fn):
        sh = kshape(Km, j['gens']); write_shape(sh, j['q'], fn)
    extra = ' '.join(j.get('extra', []))
    cmd = f'ulimit -v 3500000; ulimit -t {j["tmo"]}; exec nice -n {j.get("nice",8)} {QP} {fn} real=1 l2c=1 first=1 maxprint=1 sym=1 uA=1 tF=1 s3=1 prog=60 {extra}'
    fo = open(f'logs/{j["name"]}.out', 'w'); fe = open(f'logs/{j["name"]}.err', 'w')
    t0 = time.time()
    p = subprocess.Popen(['bash', '-c', cmd], stdout=fo, stderr=fe, preexec_fn=os.setsid, env=dict(os.environ, OMP_NUM_THREADS='1'))
    try:
        rc = p.wait(timeout=8 * j['tmo']); to = False
    except subprocess.TimeoutExpired:
        os.killpg(p.pid, signal.SIGKILL); p.wait(); rc = None; to = True
    fo.close(); fe.close()
    out = open(f'logs/{j["name"]}.out').read(); err = open(f'logs/{j["name"]}.err').read()
    sol = [l for l in out.split('\n') if l.startswith('SOL')]
    stats = [l for l in out.split('\n') if l.startswith('STATS')]
    prog = [l for l in err.split('\n') if l.startswith('PROG ')]
    lists = [l for l in err.split('\n') if l.startswith('lists:')]
    cprobes = [l for l in err.split('\n') if l.startswith('CHUNK_PROBES')]
    if (not sol) and (not stats) and rc not in (0, None): to = True
    if sol: status = 'HIT'; open(f'hits/{j["name"]}.sol', 'w').write(sol[0] + '\n')
    elif to: status = 'TIMEOUT'
    elif rc == 0 and stats:
        m = re.search(r'l2_toobig=(\d+)', stats[0]); status = 'INCOMPLETE' if (m and int(m.group(1)) > 0) else 'EMPTY'
    else: status = 'ERROR'
    rec = dict(name=j['name'], v=j['v'], K=j['K'], q=j['q'], gens=j['gens'], status=status, rc=rc, tmo=j['tmo'], sec=round(time.time() - t0, 1), extra=extra, sol=sol[0] if sol else None, stats=stats[0] if stats else None,
               lastprog=prog[-1] if prog else None, lists=lists[0] if lists else None, chunks=cprobes[:3], errtail=err.strip().split('\n')[-1][:200] if status == 'ERROR' else None)
    open(LOG, 'a').write(json.dumps(rec) + '\n')
    print(status, j['name'], rec['sec'], (stats[0][:200] if stats else (prog[-1] if prog else '')), flush=True)

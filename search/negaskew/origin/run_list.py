"""run qp over an explicit list of shapes: python3 run_list.py TMO OUT.jsonl "m:q,m:q,..." [extra qp args]   (H_K: all odd-order subgroups with nu in [NUMIN,NUMAX] given by env NUMIN/NUMAX)"""
import sys, os, json, subprocess, time, resource
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from shapes_lib import *
HERE = os.path.dirname(os.path.abspath(__file__))
TMO = int(sys.argv[1]); OUT = sys.argv[2]; pairs = [tuple(map(int, x.split(':'))) for x in sys.argv[3].split(',')]; extra = sys.argv[4:]
NUMIN = int(os.environ.get('NUMIN', 3)); NUMAX = int(os.environ.get('NUMAX', 9)); BIN = os.environ.get('QPBIN', os.path.join('..', 'qp'))      # engine, relative to this directory (qp_oo: built here) or absolute
os.makedirs(os.path.join(HERE, 'shapes'), exist_ok=True)
done = set()
if os.path.exists(OUT):
    for l in open(OUT):
        try: d = json.loads(l); done.add((d['m'], d['q'], tuple(d['hK'])))
        except Exception: pass
def lim(): resource.setrlimit(resource.RLIMIT_AS, (3500000 * 1024, 3500000 * 1024))
for m, q in pairs:
    for H, gens in sorted(subgroups_odd(m).items(), key=lambda kv: len(kv[0])):
        nu = len(orbits(m, H))
        if nu < NUMIN or nu > NUMAX: continue
        if (m, q, tuple(gens)) in done: continue
        fn = os.path.join(HERE, 'shapes', f"o_{m}_{q}_{'_'.join(map(str, gens))}.txt")
        subprocess.run([sys.executable, os.path.join(HERE, '..', 'qshape.py'), json.dumps([m]), str(q), json.dumps(gens), fn], capture_output=True, text=True, cwd=HERE)
        t0 = time.time()
        cmd = ['nice', '-n', '5', os.path.join(HERE, BIN), fn, 'l2c=1'] + extra
        try:
            p = subprocess.run(cmd, capture_output=True, text=True, timeout=TMO, preexec_fn=lim); out, err = p.stdout, p.stderr; status = 'done' if p.returncode == 0 else f'rc{p.returncode}'
        except subprocess.TimeoutExpired as e:
            out = e.stdout.decode() if e.stdout else ''; err = e.stderr.decode() if e.stderr else ''; status = 'timeout'
        cnt = None; stats = ''
        for l in out.split('\n'):
            if l.startswith('COUNT'): cnt = int(l.split()[1])
        for l in err.split('\n') + out.split('\n'):
            if l.startswith('STATS'): stats = l
        last = [l for l in (err + out).split('\n') if l.startswith('PROGRESS')][-1:]
        sols = [l for l in out.split('\n') if l.startswith('SOL')]
        rec = dict(m=m, q=q, hK=gens, nu=nu, v=m * q, status=status, count=cnt, t=round(time.time() - t0, 1), stats=stats, last=last, sols=sols[:3], extra=extra, bin=BIN)
        open(OUT, 'a').write(json.dumps(rec) + '\n')
        print(m, q, gens, 'nu', nu, status, cnt, rec['t'], stats[:140], last[:1], flush=True)

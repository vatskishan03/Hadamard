"""run qp (d1=1 origin-defect restriction) over the origin-defect shape list.  usage: python3 run_oo.py MMIN MMAX NUMIN NUMAX TMO OUTJSONL [extra qp args...]"""
import sys, os, json, subprocess, time, resource
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from oo_shapes import shapes
HERE = os.path.dirname(os.path.abspath(__file__))
mmin, mmax, NUMIN, NUMAX, TMO = (int(x) for x in sys.argv[1:6]); OUT = sys.argv[6]; extra = sys.argv[7:]
QPBIN = os.environ.get('QPBIN', os.path.join(HERE, '..', 'qp'))      # the engine: ../qp (negaskew), or an absolute path
os.makedirs(os.path.join(HERE, 'shapes'), exist_ok=True)
done = set()
if os.path.exists(OUT):
    for l in open(OUT):
        try: d = json.loads(l); done.add((d['m'], d['q'], tuple(d['hK'])))
        except Exception: pass
sh = shapes(mmin, mmax, NUMAX, NUMIN)
sh.sort(key=lambda s: (s['nu'], s['v']))
def lim():
    resource.setrlimit(resource.RLIMIT_AS, (3500000 * 1024, 3500000 * 1024))
for s in sh:
    key = (s['m'], s['q'], tuple(s['gens']))
    if key in done: continue
    fn = os.path.join(HERE, 'shapes', f"o_{s['m']}_{s['q']}_{'_'.join(map(str, s['gens']))}.txt")
    subprocess.run([sys.executable, os.path.join(HERE, '..', 'qshape.py'), json.dumps([s['m']]), str(s['q']), json.dumps(s['gens']), fn], capture_output=True, text=True, cwd=HERE)
    t0 = time.time()
    cmd = ['nice', '-n', '5', QPBIN, fn, 'l2c=1', 'd1=1', 'maxprint=3'] + extra
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=TMO, preexec_fn=lim)
        out = p.stdout; err = p.stderr; status = 'done' if p.returncode == 0 else f'rc{p.returncode}'
    except subprocess.TimeoutExpired as e:
        out = e.stdout.decode() if e.stdout else ''; err = e.stderr.decode() if e.stderr else ''; status = 'timeout'
    cnt = None; stats = ''
    for l in out.split('\n'):
        if l.startswith('COUNT'): cnt = int(l.split()[1])
    for l in err.split('\n') + out.split('\n'):
        if l.startswith('STATS'): stats = l
    last = [l for l in (err + out).split('\n') if l.startswith('PROGRESS')][-1:] 
    sols = [l for l in out.split('\n') if l.startswith('SOL')]
    rec = dict(m=s['m'], q=s['q'], hK=s['gens'], nu=s['nu'], v=s['v'], status=status, count=cnt, t=round(time.time() - t0, 1), stats=stats, last=last, sols=sols[:3], extra=extra)
    open(OUT, 'a').write(json.dumps(rec) + '\n')
    print(s['m'], s['q'], s['gens'], 'nu', s['nu'], status, cnt, rec['t'], stats[:120], flush=True)

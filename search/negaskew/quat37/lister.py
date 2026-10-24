"""qp37q listonly=1 for every identity-feasible shape of an enum37q.py file: list sizes |A|, |L| (resume-safe jsonl).  usage: python3 lister.py SHAPES.json TMO LOG.jsonl      (binary: $QP, default ./qp37q)"""
import json, subprocess, sys, os, time, re
QP = os.environ.get('QP', './qp37q')
sh=json.load(open(sys.argv[1])); TMO=int(sys.argv[2]); LOG=sys.argv[3]
done=set()
if os.path.exists(LOG):
    for l in open(LOG): d=json.loads(l); done.add((d['v'],d['K'],d['q'],tuple(d['gens'])))
os.makedirs('shapes',exist_ok=True)
for s in sorted([x for x in sh if x['ident']],key=lambda x:(x['nu'],x['v'])):
    key=(s['v'],s['K'],s['q'],tuple(s['gens']))
    if key in done: continue
    fn=f"shapes/c_{s['v']}_{s['K']}_{s['q']}_{'_'.join(map(str,s['gens']))}.txt"
    if not os.path.exists(fn):
        subprocess.run([sys.executable,'qshape37q.py',json.dumps([s['K']]),str(s['q']),json.dumps(s['gens']),fn],capture_output=True)
    t=time.time()
    try:
        r=subprocess.run(['bash','-c',f'ulimit -v 3500000; export OMP_NUM_THREADS=1; nice -n 5 {QP} {fn} listonly=1'],capture_output=True,text=True,timeout=TMO)
        m=re.search(r'LISTS \|A\|=(\d+) classesA=(\d+) \|L\|=(\d+) classesL=(\d+)',r.stdout)
        res=dict(A=int(m.group(1)),cA=int(m.group(2)),L=int(m.group(3)),cL=int(m.group(4))) if m else dict(err=(r.stderr[-200:]))
    except subprocess.TimeoutExpired:
        res=dict(err='TIMEOUT')
    d=dict(v=s['v'],K=s['K'],q=s['q'],gens=s['gens'],nu=s['nu'],ells=s['ells'],sec=round(time.time()-t,1),**res)
    open(LOG,'a').write(json.dumps(d)+'\n'); print(d,flush=True)

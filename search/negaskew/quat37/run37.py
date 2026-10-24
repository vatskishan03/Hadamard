"""exhaustive (or capped) quaternary plug-in runs. usage: python3 run37.py LISTS.jsonl OUT.jsonl TMO [nu list e.g. 5] [extra qp37q args...]   (binary: $QP, default ./qp37q)
LISTS.jsonl rows: v,K,q,gens,nu,A,L (from lister.py).  Resume-safe: skips shapes already in OUT.  Skips v with 8v-1 prime power (classical Paley I) unless INCLUDE_PALEY=1."""
import json, subprocess, sys, os, time, re
from classical import is_pp
LISTS=sys.argv[1]; OUT=sys.argv[2]; TMO=int(sys.argv[3]); NUS=[int(x) for x in sys.argv[4].split(',')] if len(sys.argv)>4 else [5]
EXTRA=sys.argv[5:]
QPBIN=os.environ.get('QP','./qp37q')
inc=os.environ.get('INCLUDE_PALEY')=='1'
done=set()
if os.path.exists(OUT):
    for l in open(OUT): d=json.loads(l); done.add((d['v'],d['K'],d['q'],tuple(d['gens'])))
rows=[json.loads(l) for l in open(LISTS)]
rows=[r for r in rows if r['nu'] in NUS and 'A' in r and r['A']>0 and r['L']>0 and (inc or not is_pp(8*r['v']-1))]
rows.sort(key=lambda r:r['A']*r['L'])
os.makedirs('runs',exist_ok=True)
for r in rows:
    key=(r['v'],r['K'],r['q'],tuple(r['gens']))
    if key in done: continue
    tag=f"{r['v']}_{r['K']}_{r['q']}_{'_'.join(map(str,r['gens']))}"
    fn=f"shapes/c_{tag}.txt"
    if not os.path.exists(fn): subprocess.run([sys.executable,'qshape37q.py',json.dumps([r['K']]),str(r['q']),json.dumps(r['gens']),fn],capture_output=True)
    t=time.time(); st='?'
    cmd=f'ulimit -v 3500000; export OMP_NUM_THREADS=1; nice -n 5 {QPBIN} {fn} l2c=1 first=1 maxprint=3 '+' '.join(EXTRA)
    try:
        p=subprocess.run(['bash','-c',cmd],capture_output=True,text=True,timeout=TMO)
        out,err=p.stdout,p.stderr
        sols=[l for l in out.splitlines() if l.startswith('SOL')]
        m=re.search(r'l1_instances=(\d+).*l2_solutions=(\d+) l2_toobig=(\d+)',out)
        if p.returncode!=0: st='ERROR rc=%d %s'%(p.returncode,err[-150:])
        elif sols: st='HIT'
        elif m and int(m.group(3))==0: st='EMPTY'
        elif m: st='INCOMPLETE_TOOBIG'
        else: st='NOSTATS'
        res=dict(l1=int(m.group(1)) if m else None,l2sol=int(m.group(2)) if m else None,toobig=int(m.group(3)) if m else None,sols=sols)
    except subprocess.TimeoutExpired as e:
        st='TIMEOUT'; res=dict(err=(e.stderr.decode()[-300:] if e.stderr else ''))
    d=dict(v=r['v'],K=r['K'],q=r['q'],gens=r['gens'],nu=r['nu'],A=r['A'],L=r['L'],status=st,sec=round(time.time()-t,1),**res)
    open(OUT,'a').write(json.dumps(d)+'\n'); print(json.dumps(d)[:300],flush=True)

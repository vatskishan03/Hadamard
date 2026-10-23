"""Runs the spin engine on the shapes of a list written by enum_spin.py (or spin_filter.py) with NLO <= n <= NHI orbits;
resumable (skips the shapes already in the output file).
usage: python3 spin_sweep.py NLO NHI seconds out.jsonl [shapes.json]      (shapes.json default runs/spin_shapes_n17.json)
Writes one json line per shape: status HIT / empty (exhaustive) / timeout / error, seconds, and the SOL or COUNT line.
For n >= 15 first remove the shapes that fail the row-sum test with spin_filter.py (runs/spin_rank_dc.json): a shape
without a row-sum pattern is empty at once.  Example (NOTES.md section 88), after enum_spin.py 17, spin_est.py and spin_filter.py:
python3 spin_sweep.py 14 17 3000 runs/spin_sweep17.jsonl runs/spin_rank_dc.json    both shapes (v = 817, H = <814> and <769>, mu = 4) are
EXHAUSTIVELY empty: COUNT 0 (passDC 6830208)."""
import json,sys,time,os,subprocess
from spin_drv import *
nlo=int(sys.argv[1]); nhi=int(sys.argv[2]); tmo=int(sys.argv[3]); outf=sys.argv[4]
shapes=json.load(open(sys.argv[5] if len(sys.argv)>5 else 'runs/spin_shapes_n17.json'))
done=set()
if os.path.exists(outf):
    for l in open(outf):
        r=json.loads(l); done.add((r['v'],tuple(r['gens']),r['mu']))
for row in shapes:
    n,v,m,o,g,mu,K=row[-7:]      # rows of enum_spin.py have 7 entries, rows of spin_est.py / spin_filter.py have the estimate in front
    if not(nlo<=n<=nhi) or (v,tuple(g),mu) in done: continue
    t0=time.time()
    try:
        d,op,out,err=run_spin(v,g,mu,maxprint=1,timeout=tmo)
        sol=[l for l in out.split('\n') if l.startswith('SOL')]
        cnt=[l for l in out.split('\n') if l.startswith('COUNT')]
        status='HIT' if sol else ('empty' if cnt else 'error')
    except subprocess.TimeoutExpired:
        sol=[];cnt=[];status='timeout'
    rec=dict(v=v,m5=m,n=n,Hord=o,gens=g,mu=mu,Kord=K,status=status,sec=round(time.time()-t0,1),sol=sol[0] if sol else None,cnt=cnt[0] if cnt else None)
    os.makedirs(os.path.dirname(outf) or '.',exist_ok=True)
    open(outf,'a').write(json.dumps(rec)+'\n'); print(rec,flush=True)

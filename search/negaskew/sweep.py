import json,sys,time,os,subprocess
from qdrv import *
nmax=int(sys.argv[1]); vmax=int(sys.argv[2]); tmo=int(sys.argv[3]); outf=sys.argv[4]; mode=int(sys.argv[5]) if len(sys.argv)>5 else 1
cands=json.load(open(f'enum_n{nmax}_v{vmax}.json')) if os.path.exists(f'enum_n{nmax}_v{vmax}.json') else json.load(open('enum_n10_v1500.json'))
done=set()
if os.path.exists(outf):
    for l in open(outf):
        r=json.loads(l); done.add((r['v'],tuple(r['gens']),r['mode']))
jobs=[]
for v,(m,sh) in cands.items():
    v=int(v)
    if v>vmax: continue
    for n,o,g in sh:
        if n<=int(os.environ.get("NRUN",nmax)): jobs.append((n,v,m,o,g))
jobs.sort()
for n,v,m,o,g in jobs:
    if (v,tuple(g),mode) in done: continue
    t0=time.time()
    fn=f"shapes/sw_{v}_{n}_{o}.txt"; d=shape(v,g,fn)
    try:
        r=subprocess.run(['./qskew2',fn,str(mode),'1','1'],capture_output=True,text=True,timeout=tmo)
        out=r.stdout; status='exhaustive_empty'
        hit=[l for l in out.split('\n') if l.startswith('SOL')]
        if hit: status='HIT'
        elif 'COUNT' not in out: status='error'
        else:
            inc=[int(l.split()[1]) for l in out.split('\n') if l.startswith('INCOMPLETE_TUPLES')]
            if inc and inc[0]>0: status='incomplete_%d'%inc[0]
    except subprocess.TimeoutExpired as ex:
        out=''; hit=[]; status='timeout'
    rec=dict(v=v,m5=m,n=n,order=o,gens=g,mode=mode,status=status,sec=round(time.time()-t0,1),sol=hit[0] if hit else None)
    open(outf,'a').write(json.dumps(rec)+'\n')
    print(rec,flush=True)

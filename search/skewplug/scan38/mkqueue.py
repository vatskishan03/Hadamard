"""per scanned shape: the best slots by estimated completions per class (logrho), the cost of a slot in classes N1 / (2 (4 |Fix|)^3) and the sums over the slots; -> runs/shapes_scanned.json.  Reads runs/table38.json (table38.py) and est/S*.ext*.jsonl (scan38.py).  usage (in search/skewplug): python3 scan38/mkqueue.py"""
import json, glob, math, os, collections
os.makedirs('runs', exist_ok=True)
tab=json.load(open('runs/table38.json'))
print(collections.Counter(t['status'] for t in tab))
out=[]
for t in tab:
    if t['status']!='SCANNED': continue
    files=sorted(glob.glob(f"est/S{t['idx']:04d}_*.ext*.jsonl"))
    if not files: continue
    for f in files:
        lo=f.endswith('_lo.ext.jsonl')
        rows=[]
        for l in open(f):
            try: r=json.loads(l)
            except: continue
            if 'logE' in r and r['logE']>-50 and 'nopat' not in r: rows.append(r)
        if not rows: continue
        N=1
        for m in t['kmods']: N*=m
        # fixed points of H_K (cyclic only)
        cyc=len(t['kmods'])==1
        F=None
        if cyc:
            n=t['kmods'][0]; hs=[g[0] for g in t['gens']]
            F=sum(1 for k in range(n) if all((h-1)*k%n==0 for h in hs))
        for r in rows: r['cost']=10**r['logN1']/(2*(4*(F or 1))**3) if cyc else None
        rows.sort(key=lambda r:-r['logrho'])
        top=[r for r in rows if r['logE']>-4][:30]
        best=top[0] if top else rows[0]
        out.append(dict(idx=t['idx'],v=t['v'],q=t['q'],kmods=t['kmods'],gens=t['gens'],order=t['order'],nu=t['nu'],lo=lo,nfeas=len(rows),F=F,file=os.path.basename(f),
                        best_rho=best['logrho'],best_rho_E=best['logE'],best_cost=best['cost'],
                        cum_cost30=sum(r['cost'] for r in top if r['cost']) if cyc else None, sumE=sum(10**r['logE'] for r in rows)))
json.dump(out,open('runs/shapes_scanned.json','w'))
out=[o for o in out if not o['lo']]
out.sort(key=lambda o:-o['best_rho'])
print(len(out))
for o in out[:60]: print(o['idx'],o['v'],o['q'],o['kmods'],o['gens'],'o',o['order'],'nu',o['nu'],'nfeas',o['nfeas'],'F',o['F'],'best rho %.2f (logE %.2f) cost %.3g'%(o['best_rho'],o['best_rho_E'],o['best_cost'] or -1),'cum30 %.3g'%(o['cum_cost30'] or -1),'sumE %.3g'%o['sumE'])

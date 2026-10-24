"""the shapes of runs/shapes_scanned.json ranked by the summed estimate E of their slots of cost <= CAP classes (env CAP, default 5000).  usage (in search/skewplug): python3 scan38/mkqueue2.py"""
import json, glob, math, os
sc=json.load(open('runs/shapes_scanned.json'))
CAP=float(os.environ.get('CAP','5000'))
rows=[]
for o in sc:
    if o['lo'] or o['F'] is None: continue
    recs=[]
    for l in open('est/'+o['file']):
        try: r=json.loads(l)
        except: continue
        if 'logE' in r and r['logE']>-50 and 'nopat' not in r: recs.append(r)
    F=o['F']
    cheap=[r for r in recs if 10**r['logN1']/(2*(4*F)**3) <= CAP]
    cost=sum(10**r['logN1']/(2*(4*F)**3) for r in cheap); E=sum(10**r['logE'] for r in cheap)
    allE=sum(10**r['logE'] for r in recs)
    rows.append((E,o['idx'],o['v'],o['q'],o['kmods'],o['gens'],o['order'],len(recs),len(cheap),cost,allE))
rows.sort(reverse=True)
for r in rows[:50]: print('cheapE %.3g idx %d v %d q %d K %s gens %s o %d nfeas %d ncheap %d cheapcost %.3g allE %.3g'%r)

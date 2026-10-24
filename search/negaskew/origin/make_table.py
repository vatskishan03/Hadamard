import json, os
rows=[]
def l1(st):
    return int(st.split('l1_instances=')[1].split()[0]) if 'l1_instances=' in st else None
for fn,tag in (('data/runs_oo5.jsonl','oo5'),('data/runs_oo7.jsonl','oo7'),('data/runs_oo7long.jsonl','oo7long'),('data/runs_small.jsonl','small')):
    if not os.path.exists(fn): continue
    for l in open(fn):
        d=json.loads(l)
        verdict={'done':('HIT' if (d['count'] or 0)>0 else 'EMPTY (exhaustive in class)'),'timeout':'UNKNOWN (timeout)'}.get(d['status'],d['status'])
        if d['status']=='done' and (d['count'] or 0)>0 and 'first=1' in ' '.join(d.get('extra',[])): verdict='HIT (first=1 stop)'
        rows.append((d['nu'],d['v'],d['m'],d['q'],str(d['hK']),verdict,l1(d.get('stats','')),d['t'],tag))
rows.sort()
with open('runs/results_table.txt','w') as f:
    f.write('nu  v   m   q  H_K-gens  verdict  level1_instances  seconds  run\n')
    for r in rows: f.write('%d %d %d %d %s %s %s %.0f %s\n'%r)
print(len(rows),'rows')

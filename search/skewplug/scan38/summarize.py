"""per exh38.py job (runs/*.slots.jsonl): slots, finished slots, slots skipped (A list too large, join cap), classes sent to level 2, completions, hits, CPU seconds.  usage (in search/skewplug): python3 scan38/summarize.py"""
import json, glob, os, collections, re
W='runs/'
rows=[]
for f in sorted(glob.glob(W+'*.slots.jsonl')):
    name=os.path.basename(f)[:-len('.slots.jsonl')]
    if name.startswith('tst') or name.startswith('hit'): continue
    recs=[]
    for l in open(f):
        try: recs.append(json.loads(l))
        except: pass
    last={}
    for r in recs: last[(tuple(r['al']),tuple(r['ms']))]=r      # last record per slot
    fin=sum(1 for r in last.values() if r.get('finished'))
    nskipA=sum(1 for r in last.values() if r.get('skip')=='A_toobig'); nskipJ=sum(r.get('nskip',0) for r in last.values())
    inc=sum(1 for r in last.values() if not r.get('complete',True))
    nn=sum(r.get('nnew',0) for r in recs); nk=sum(r.get('nok',0) for r in recs); t=sum(r.get('t_slot',0) for r in recs)
    hits=sum(len(r.get('HIT',[])) for r in recs)
    seen=os.path.getsize(W+name+'.seen')//8 if os.path.exists(W+name+'.seen') else 0
    rows.append((name,len(last),fin,nskipA,nskipJ,inc,nn,nk,hits,round(t),seen))
print('%-14s %6s %6s %6s %6s %6s %9s %5s %5s %7s %8s'%('job','slots','finish','A_big','jskip','nocomp','classes','compl','hits','cpu_s','seenkeys'))
for r in rows: print('%-14s %6d %6d %6d %6d %6d %9d %5d %5d %7d %8d'%r)

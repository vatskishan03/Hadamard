"""low-energy slot files: all (allocation, m-split) slots of E < (lowest energy of the scan's est file), energy ascending; est fields are dummies (logN1 = 0).  usage (in search/skewplug): python3 scan38/mklow.py IDX ...  -> est/low<idx>_v<v>.jsonl   (IDX of scan38/data/queue/shapes.json; needs the slot file s['est'] of the shape, a copy of its scan file est/S*.ext.jsonl)"""
import sys, json, os, math, random
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..'))
from sk import SShape
from sk_pat import allocations, m_splits
shapes = {s['idx']: s for s in json.load(open(os.environ.get('ROT_QUEUE', os.path.join(HERE, 'data', 'queue', 'shapes.json'))))}      # env ROT_QUEUE: another shape list (as for rotate.py)
for idx in [int(x) for x in sys.argv[1:]]:
    s = shapes[idx]; sh = SShape(s['K'], s['q'], s['hu'])
    est = [json.loads(l) for l in open(s['est'])]
    emin_est = min(sum(r['al']) for r in est)
    al = allocations(sh); slots = []
    for a, em in al:
        E = sum(a)
        if E >= emin_est: continue
        for ms in m_splits(a, em): slots.append((E, list(a), list(ms)))
    slots.sort(key=lambda x: (x[0], x[1], x[2]))
    os.makedirs('est', exist_ok=True); path = f"est/low{idx}_v{s['v']}.jsonl"
    with open(path, 'w') as f:
        for E, a, ms in slots: f.write(json.dumps(dict(al=a, ms=ms, logN1=0.0, logE=0.0, logrho=0.0)) + '\n')
    print(idx, s['v'], 'E <', emin_est, 'slots', len(slots), path)

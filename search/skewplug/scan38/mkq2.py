"""queue of the deeper scan passes: the shapes with margin >= MN (default -4) whose scan is not COMPLETE, best first -> runs/queue2.json.  usage (in search/skewplug): python3 scan38/mkq2.py [MN]"""
import json, sys, os
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data'); os.makedirs('runs', exist_ok=True)
mn = float(sys.argv[1]) if len(sys.argv) > 1 else -4.0
S = json.load(open('runs/summary38.json')); shapes = json.load(open(DATA + '/shapes38.json'))
tab = {}
for s in S:
    if s.get('window', 'top') != 'top': continue
    tab[s['idx']] = s
sel = []
for i, s in tab.items():
    r = shapes[i]
    if r['decided']: continue
    cov = s.get('coverage')
    if s.get('margin') is not None and s['margin'] >= mn and s['status'] != 'COMPLETE': sel.append((s['margin'], i))
sel.sort(reverse=True)
q = [i for _, i in sel]
json.dump(q, open('runs/queue2.json', 'w')); print(len(q), q)

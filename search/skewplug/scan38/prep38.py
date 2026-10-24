"""annotate the output of enum38.py: decided-match (exhaustive verdicts of NOTES 88-90), identity feasibility (allocations), pool size.  usage (in search/skewplug): python3 scan38/prep38.py RAW OUT      (reads runs/decided38.json)"""
import sys, os, json, time, warnings
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..'))      # scan38/ and the modules of search/skewplug
DATA = os.path.join(HERE, 'data')
from pshape import shape_from_row
from sk_pat import allocations, m_splits
import numpy as np
from awin import alloc_window
dec = json.load(open('runs/decided38.json'))
dkey = {}; dnc = {}
for d in dec:
    if d['kind'] == 'cyc':
        dkey[(d['v'], d['q'], d['m'], frozenset(tuple(g) for g in d['group']))] = d['src']
    dnc.setdefault((d['v'], d['q'], d['m'], d['order'], d['nu'], tuple(d['sizes'])), d['src'])
rows = json.load(open(sys.argv[1]))
out = []
t0 = time.time()
for i, r in enumerate(rows):
    r['idx'] = i
    blocks = r['blocks']
    r['decided'] = None
    if len(blocks) == 1 and blocks[0][0] == 'cyc':
        m = blocks[0][1]
        G = frozenset(tuple(g[0] for g in [()]) if False else (g,) for g in []) 
        # group elements from gens: closure over units
        gens = [g[0] for g in r['gens']]
        S = {1}; fr = [1]
        while fr:
            x = fr.pop()
            for g in gens:
                y = x * g % m
                if y not in S: S.add(y); fr.append(y)
        key = (r['v'], r['q'], m, frozenset((u,) for u in S))
        if key in dkey: r['decided'] = dkey[key]
    else:
        key = (r['v'], r['q'], r['m'], r['order'], r['nu'], tuple(r['sizes']))
        if key in dnc: r['decided'] = 'probable: ' + dnc[key]
    try:
        sh = shape_from_row(r)
        win, emin, emax = alloc_window(sh, 5)
        r['nwin_alloc'] = len(win); r['Emin'] = emin; r['Emax'] = emax
        r['nalloc'] = len(win)          # (window allocations; 0 <=> the identity cannot be realised by the orbit weights at all)
        tot = 0
        for (a, em) in win:
            poly = np.array([1.0])
            for j, n in enumerate(a):
                p = np.zeros(n + 1)
                for mm in range(0, n + 1, 2 if j == 0 else 1): p[mm] = 1
                poly = np.convolve(poly, p)
            tot += float(poly[em]) if em < len(poly) else 0.0
        r['pool5_ub'] = int(tot) if win else 0
    except Exception as ex:
        r['err'] = str(ex)[:100]
    out.append(r)
    if i % 25 == 0: print(i, round(time.time() - t0), flush=True)
json.dump(out, open(sys.argv[2], 'w'))
from collections import Counter
print('decided', sum(1 for r in out if r['decided']), 'dead identity', sum(1 for r in out if r.get('nalloc') == 0), 'err', sum(1 for r in out if 'err' in r))

"""re-annotate the decided field of an enumeration (enum38.py output) after decided38.py was re-run: python3 scan38/flag38.py IN.json OUT.json"""
import json, sys
dec = json.load(open('runs/decided38.json'))
dkey = {}; dnc = {}
for d in dec:
    if d['kind'] == 'cyc':
        dkey[(d['v'], d['q'], d['m'], frozenset(tuple(g) for g in d['group']))] = d['src']
    dnc.setdefault((d['v'], d['q'], d['m'], d['order'], d['nu'], tuple(d['sizes'])), d['src'])
rows = json.load(open(sys.argv[1]))
for i, r in enumerate(rows):
    r['idx'] = i; r['decided'] = None
    if len(r['blocks']) == 1 and r['blocks'][0][0] == 'cyc':
        m = r['blocks'][0][1]; gens = [g[0] for g in r['gens']]
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
json.dump(rows, open(sys.argv[2], 'w'))
from collections import Counter
und = [r for r in rows if not r['decided']]
print('decided', len(rows) - len(und), 'undecided', len(und), sorted(Counter(r['nu'] for r in und).items()))
print('decided by nu', sorted(Counter(r['nu'] for r in rows if r['decided']).items()))

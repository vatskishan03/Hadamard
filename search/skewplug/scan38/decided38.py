"""collect the shapes already decided EXHAUSTIVELY (EMPTY / EXACTLY EMPTY / HIT) in the sweeps of search/NOTES.md sections 88-90.  decided38.json: list of dicts
   dict(v, q, m, kind='cyc'|'nc', units=[...gens], mods=[...], src, status, order, nu, sizes)
usage (in search/skewplug): python3 scan38/decided38.py      -> runs/decided38.json   (reads scan38/data: census_shapes15.json, timeouts.json, SHAPE_TABLE.txt, verdict_table_cyclic.txt and ../negaskew/real37/data/RESULTS.txt)"""
import re, json, os, sys, itertools
from math import gcd
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..'))      # scan38/ and the modules of search/skewplug
DATA = os.path.join(HERE, 'data')
from pk import KGrp

def orbit_data(mods, units):
    Kg = KGrp(mods); N = Kg.n
    perms = [Kg.mulperm(tuple(u)) for u in units]
    orb = -np.ones(N, dtype=np.int64); n = 0
    for k0 in range(N):
        if orb[k0] >= 0: continue
        orb[k0] = n; st = [k0]
        while st:
            k = st.pop()
            for P in perms:
                k2 = int(P[k])
                if orb[k2] < 0: orb[k2] = n; st.append(k2)
        n += 1
    sizes = sorted(np.bincount(orb).tolist())
    # group order = closure size
    e = tuple([1] * len(mods)); G = {e}; fr = [e]
    while fr:
        x = fr.pop()
        for u in units:
            y = tuple((a * b) % mm for a, b, mm in zip(x, u, mods))
            if y not in G: G.add(y); fr.append(y)
    return n, sizes, len(G), frozenset(G)

out = []
def add(v, q, mods, units, src, status):
    m = int(np.prod(mods))
    if units and not isinstance(units[0], (list, tuple)): units = [[u] * len(mods) for u in units]
    nu, sizes, order, G = orbit_data(mods, units)
    out.append(dict(v=v, q=q, m=m, kind='cyc' if len(mods) == 1 else 'nc', mods=list(mods), units=[list(u) for u in units], group=sorted(map(list, G)), src=src, status=status, order=order, nu=nu, sizes=sizes))

# ---- r36 census sweep (real class, K = Z_m): every shape of census_shapes15.json except the r36 timeouts
shapes = json.load(open(DATA + '/census_shapes15.json'))
to = json.load(open(DATA + '/timeouts.json'))
tkeys = {(t['v'], t['K'], t['q'], tuple(t['gens'])) for t in to}
for (nu, v, m, q, gens, hs) in shapes:
    if (v, m, q, tuple(gens)) in tkeys: continue
    add(v, q, [m], [[g] for g in gens], 'r36 real sweep (exhaustive within time limits, not a timeout)', 'EMPTY')

# ---- r37 RESULTS (real class)
parts = {}
for l in open(os.path.join(HERE, '..', '..', 'negaskew', 'real37', 'data', 'RESULTS.txt')):
    mt = re.match(r'(EMPTY|HIT|TIMEOUT)\s+v=(\d+)\s+(\S+)\s+q=(\d+)', l)
    if not mt: continue
    st, v, name, q = mt.group(1), int(mt.group(2)), mt.group(3), int(mt.group(4))
    if st != 'EMPTY': continue
    nm = re.match(r'[cn](\d+)_K(\S+?)_q(\d+)_h(.*?)(_p\d)?$', name)
    if not nm: print('unparsed', name, file=sys.stderr); continue
    Ks, hs = nm.group(2), nm.group(4)
    mods = [int(x) for x in Ks.split('x')]
    if hs.startswith('u'):
        units = [[int(x) for x in hs[1:].split('.')]]
    else:
        units = [[int(x)] * len(mods) for x in hs.split('_')]
    key = (v, q, tuple(mods), json.dumps(units))
    parts.setdefault(key, []).append(nm.group(5))
for (v, q, mods, ju), ps in parts.items():
    if any(p for p in ps) and len([p for p in ps if p]) < 2: continue      # only one half of a split job
    add(v, q, list(mods), json.loads(ju), 'r37 real exhaustive (qp37 with symmetry reductions)', 'EMPTY')

# ---- r35 SHAPE_TABLE (statuses EXACTLY EMPTY / EXHAUSTIVE...)
for l in open(DATA + '/SHAPE_TABLE.txt'):
    mt = re.match(r'v=\s*(\d+) q=\s*(\d+) K=(\S+)\s+H_K=(\[\[.*?\]\])\s+ord=\s*(\d+) nu=\s*(\d+)', l)
    if not mt: continue
    v, q = int(mt.group(1)), int(mt.group(2))
    mods = [int(x) for x in mt.group(3).split('x')]
    rest = l[mt.end():]
    if not re.search(r'EXACTLY EMPTY|EXHAUSTIVELY EMPTY|EXHAUSTIVE per pattern \(.*?0 completions', rest): continue
    if 'UNKNOWN' in rest: continue
    H = json.loads(mt.group(4).replace(' ', ''))
    if len(mods) == 1: units = [[h[0]] for h in H]
    else: units = [list(h) for h in H]
    try: add(v, q, mods, units, 'r35 SHAPE_TABLE: ' + rest.strip()[:80], 'EMPTY')
    except Exception as e: print('r35 fail', l, e, file=sys.stderr)

# ---- r34 verdict table (cyclic)
for l in open(DATA + '/verdict_table_cyclic.txt'):
    mt = re.match(r'v=(\d+) q=(\d+) K=Z(\d+) H_K=<(\d+)> order (\d+) nu=(\d+) : (.*)', l)
    if not mt: continue
    if 'EMPTY' not in mt.group(7): continue
    add(int(mt.group(1)), int(mt.group(2)), [int(mt.group(3))], [[int(mt.group(4))]], 'r34 sweep: ' + mt.group(7)[:60], 'EMPTY')

os.makedirs('runs', exist_ok=True); json.dump(out, open('runs/decided38.json', 'w'))
from collections import Counter
print(len(out), 'decided shape records', Counter(o['src'].split(' ')[0] for o in out), Counter(o['kind'] for o in out))

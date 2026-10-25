"""Tied census for the complex Yamada-Kiyasu targets. Restart-safe log (one line per case).
Usage: census.py maxC maxAB log [spec:b ...]   (b = 0 plain, 1 bordered)"""
import sys, os, time, json, numpy as np
from gengroup import *
maxC, maxAB, log = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
specs = [s.split(':') for s in sys.argv[4:]]
done = set()
if os.path.exists(log):
    for line in open(log):
        if line.startswith('CASE '): done.add(line.split(' | ')[0])
out = open(log, 'a')
for spec, b in specs:
    bordered = b == '1'
    t0 = time.time()
    G = Group(spec)
    subs = G.subgroups()
    # coarsening menus: quadratic / quartic classes on every field factor of odd order
    def gens_for(k):
        res = []
        for i, f in enumerate(G.f):
            if f.kind == 'F' and (f.q - 1) % k == 0:
                e = list(G.mident()); e[i] = k; res.append(tuple(e))
        return res
    menus = [[], gens_for(2), gens_for(4)]
    menus = [x for i, x in enumerate(menus) if x not in menus[:i]]
    ncase = 0
    for H, gens in sorted(subs.items(), key=lambda kv: -len(kv[0])):
        sigs = [G.perm(h) for h in gens]
        OC, _ = G.orbits(sigs)
        if len(OC) > maxC or len(OC) < 4: continue
        for cA in menus:
            for cB in menus:
                key = 'CASE %s b=%d H=%s cA=%s cB=%s' % (spec, bordered, gens, cA, cB)
                if key in done: continue
                t1 = time.time()
                try:
                    res, bits = tied_run(G, gens, bordered, cA, cB, maxC, maxAB, stop_at=3)
                except MemoryError:
                    res, bits = 'MEM', None
                line = '%s | |H|=%d bits=%s res=%s t=%.0f' % (key, len(H), bits, None if res is None or res == 'MEM' else res[0], time.time() - t1)
                if res not in (None, 'MEM') and res[0]:
                    line += ' FIRST ' + json.dumps({'n': res[1][0], 'check': res[1][1], 'seqs': res[1][2]})
                    print('FOUND', key, bits, res[0], res[1][1], flush=True)
                out.write(line + '\n'); out.flush(); ncase += 1
    out.write('SPECDONE %s b=%d cases=%d t=%.0f\n' % (spec, bordered, ncase, time.time() - t0)); out.flush()
    print('specdone', spec, bordered, ncase, '%.0f s' % (time.time() - t0), flush=True)
print('census done', flush=True)

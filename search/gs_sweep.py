"""Exhaustive Goethals-Seidel sweep over small multiplier structures of all open orders."""
import re, sys, time, json, os
from math import gcd
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gs_exhaustive as E

def subgroups(n):
    units = [u for u in range(1, n) if gcd(u, n) == 1]
    def closure(gens):
        S, todo = {1}, [1]
        while todo:
            x = todo.pop()
            for g in gens:
                y = x * g % n
                if y not in S: S.add(y); todo.append(y)
        return frozenset(S)
    cyc = {}
    for u in units: cyc.setdefault(closure([u]), u)
    subs = {H: (g,) for H, g in cyc.items()}
    items = list(cyc.items())
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            subs.setdefault(closure([items[i][1], items[j][1]]), (items[i][1], items[j][1]))
    return subs

fr = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'FRONTIER.md')).read()
ns = [int(n) for _, n in re.findall(r'^\|\s*(\d+)\s*\|\s*(\d+)\s*\|', fr, re.M)]
MIN_O, MAX_O, CAP = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
only = set(map(int, sys.argv[4:]))
log = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'gs_sweep_hits.json'), 'a')
for n in ns:
    if only and n not in only: continue
    t0 = time.time(); tried = 0
    for H, gens in sorted(subgroups(n).items(), key=lambda kv: -len(kv[0])):
        if (n - 1) in H: continue                         # -1 in H: symmetric blocks, skip
        seen = set(); O = 0
        for t in range(n):
            if t not in seen:
                seen |= {t * h % n for h in H}; O += 1
        if not (MIN_O <= O <= MAX_O): continue
        tried += 1
        found, st = E.run(n, list(gens), max_pairs=CAP, K=3, verbose=False)
        if found:
            for q, seqs, had in found:
                print(f"*** n={n} gens={gens} pattern={q}: H({4*n}) Hadamard={had}", flush=True)
                log.write(json.dumps(dict(n=n, gens=list(gens), pattern=list(q), hadamard=bool(had),
                    seqs=[''.join('+' if x > 0 else '-' for x in s) for s in seqs])) + "\n"); log.flush()
    print(f"n={n}: {tried} structures ({time.time()-t0:.0f}s)", flush=True)

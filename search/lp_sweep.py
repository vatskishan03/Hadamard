"""Legendre-pair sweep over all open orders: l = 2n - 1, every multiplier subgroup with
MIN_O <= orbits <= MAX_O, exhaustive join.  Prints every structure tried and any verified pair."""
import os
import re, sys, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lp_search as LP

fr = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'FRONTIER.md')).read()
ns = [int(n) for _, n in re.findall(r'^\|\s*(\d+)\s*\|\s*(\d+)\s*\|', fr, re.M)]
MIN_O, MAX_O = int(sys.argv[1]), int(sys.argv[2])
only = set(map(int, sys.argv[3:]))
hits = []
for n in ns:
    if only and n not in only: continue
    l = 2 * n - 1
    t0 = time.time(); tried = 0
    for H, gens in sorted(LP.subgroups(l).items(), key=lambda kv: -len(kv[0])):
        orb, sizes, S, reps = None, None, None, None
        # cheap orbit count
        seen = set(); O = 0
        for t in range(l):
            if t not in seen:
                seen |= {t * h % l for h in H}; O += 1
        if not (MIN_O <= O <= MAX_O): continue
        tried += 1
        sols, O2, N = LP.find_pairs(l, H, maxpool=3e7)
        if sols:
            a, b = LP.to_sequences(l, H, sols[0])
            ok = LP.verify_lp(a, b)
            print(f"*** n={n} l={l} |H|={len(H)} gens={gens} orbits={O}: Legendre pair found, verified={ok}", flush=True)
            hits.append(dict(n=n, l=l, gens=list(gens), a=''.join('+' if x > 0 else '-' for x in a),
                             b=''.join('+' if x > 0 else '-' for x in b), verified=bool(ok)))
            json.dump(hits, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lp_hits.json'), 'w'))
    print(f"n={n} l={l}: {tried} structures searched ({time.time()-t0:.0f}s)", flush=True)
print("HITS:", len(hits))

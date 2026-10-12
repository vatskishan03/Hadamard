"""Plain and twisted multiplier-structure Legendre pair search for given lengths l (exhaustive per structure)."""
import sys, time, json
import numpy as np
sys.path.insert(0, '/home/user/Hadamard/attack')
import lp_search as LP, tlp_model as TM, tlp_search as TS
maxo = int(sys.argv[1]); ls = [int(a) for a in sys.argv[2:]]
for l in ls:
    t0 = time.time(); tried = 0
    for H, gens in sorted(LP.subgroups(l).items(), key=lambda kv: -len(kv[0])):
        lab, k = TM.orbit_labels(l, H)
        if k > maxo: continue
        sols, O, N = LP.find_pairs(l, H, maxpool=2e8, limit=3)
        if sols is None: continue
        tried += 1
        if sols:
            a, b = LP.to_sequences(l, H, sols[0]); ok = LP.verify_lp(a, b)
            print(f"*** PLAIN l={l} |H|={len(H)} gens={gens} orbits={k}: LP verified={ok}", flush=True)
            with open('lp_open_hits.json', 'a') as f:
                f.write(json.dumps(dict(l=l, kind='plain', gens=list(gens), a=''.join('+' if x > 0 else '-' for x in a), b=''.join('+' if x > 0 else '-' for x in b), verified=bool(ok))) + "\n")
    print(f"l={l}: plain structures {tried} ({time.time()-t0:.0f}s)", flush=True)
    t0 = time.time(); tried = 0
    for H, gens, mu, lab, k in TM.twisted_structures(l, maxo):
        hits, total = TS.search(l, H, mu, lab, k, 2e8)
        if hits is None: continue
        tried += 1
        for y in hits[:1]:
            a = y[lab]; b = a[(np.arange(l) * mu) % l]; ok = TS.verify_pair(l, a, b)
            print(f"*** TWISTED l={l} |H|={len(H)} gens={gens} mu={mu} orbits={k}: LP verified={ok}", flush=True)
            with open('lp_open_hits.json', 'a') as f:
                f.write(json.dumps(dict(l=l, kind='twisted', gens=list(gens), mu=mu, a=''.join('+' if x > 0 else '-' for x in a), b=''.join('+' if x > 0 else '-' for x in b), verified=bool(ok))) + "\n")
    print(f"l={l}: twisted structures {tried} ({time.time()-t0:.0f}s)", flush=True)

"""Sweep skew Goethals-Seidel multiplier structures exactly, without trusting the Gaussian model.

For each v, each subgroup H of Z_v^* with -1 not in H and at most MAXO orbits, and each admissible row-sum
pattern, run the exact four-list matcher find4b when the pair work N_A N_D + N_B N_C is at most CAP.
Stops at the first v-structure with a verified hit (one order is enough) and records hits in skew_hits.json.

  python3 skew_sweep.py MAXO CAP v ...
"""
import json, sys, time
from pathlib import Path as _P; sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
from gs_sweep_subgroups import subgroups
import skew_exhaustive as SE
HERE = str(W.HERE)


def structures(v, maxo):
    out = []
    for H, gens in subgroups(v).items():
        if (v - 1) in H: continue
        O = len({min(t * h % v for h in H) for t in range(v)})
        if 6 <= O <= maxo: out.append((O, list(gens), len(H)))
    return sorted(out)


if __name__ == '__main__':
    maxo, cap = int(sys.argv[1]), float(sys.argv[2])
    for v in map(int, sys.argv[3:]):
        for O, gens, h in structures(v, maxo):
            t = time.time()
            try:
                found = SE.run(v, gens, cap=cap, verbose=False)
            except AssertionError:
                continue
            print(f"v={v} |H|={h} gens={gens} orbits={O}: {len(found)} verified [{time.time()-t:.0f}s]", flush=True)
            if found:
                with open(W.hits(), 'a') as f:
                    for x in found: f.write(json.dumps(x) + '\n')
                print(f"HIT v={v} 4v={4*v} gens={gens}", flush=True)
                break

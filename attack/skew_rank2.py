"""Rank skew multiplier structures for composite v with the compression filter.

E' (skew_rank.py) assumes uniform lists of at most L sequences. Restricting each sampled list to sequences
whose Z_d-compression norms occur in a feasible combination (skew_compress.role_filters) loses no
solutions, so the yield of a sampled list rises by 1/keep. Structures that fail a compression test are
dropped.

  python3 skew_rank2.py L MAXO v ...
"""
import math, sys
from math import log10
import numpy as np
sys.path.insert(0, '/home/user/Hadamard/attack')
import spence_lib as L
from gs_sweep_subgroups import subgroups
import skew_rank as R
import skew_compress as SC


def primes(v):
    return [p for p in range(2, v) if v % p == 0 and all(p % q for q in range(2, int(p ** 0.5) + 1))]


if __name__ == '__main__':
    cap = float(sys.argv[1]); maxo = int(sys.argv[2]); ns = [int(a) for a in sys.argv[3:]]
    rng = np.random.default_rng(1)
    for v in ns:
        ds = [d for d in primes(v) if d < v]
        if not ds: continue
        for H, gens in subgroups(v).items():
            if (v - 1) in H: continue
            O = len({min(t * h % v for h in H) for t in range(v)})
            if O > maxo or O < 6: continue
            st, res = R.rank(v, list(gens), cap, rng=rng)
            for lE2, lE, q, lN in sorted(res, reverse=True)[:3]:
                if lE < -2: continue
                gain, dead, info = 0.0, False, []
                for d in ds:
                    keep = SC.role_filters(st, v, d, list(q), cap=200_000)
                    if keep is None: info.append(f"d={d}:big"); continue
                    if not all(keep): dead = True; info.append(f"d={d}:infeasible"); break
                    Ys = [None] + [st.sample(r, 3000, rng) for r in q]
                    # skew A sample
                    rep = np.array([np.flatnonzero(st.orb == o)[0] for o in range(st.O)])
                    partner = st.orb[st.G.neg()[rep]]
                    pairs = [(o, partner[o]) for o in range(1, st.O) if o < partner[o]]
                    YA = np.zeros((3000, st.O)); YA[:, 0] = 1
                    sg = np.where(rng.random((3000, len(pairs))) < 0.5, 1.0, -1.0)
                    for j, (o, p_) in enumerate(pairs): YA[:, o] = sg[:, j]; YA[:, p_] = -sg[:, j]
                    Ys[0] = YA
                    fr = []
                    for role, Y in enumerate(Ys):
                        Nn = SC.norm_rows(Y, st, d)
                        fr.append(max(np.mean([tuple(x) in keep[role] for x in Nn]), 1 / 3000))
                    g = sum(-log10(f) for f, ln in zip(fr, lN) if ln > log10(cap))
                    gain += g; info.append(f"d={d}:keep=({','.join(f'{f:.3f}' for f in fr)})")
                if dead:
                    print(f"   dead n={v} |H|={len(H)} gens={','.join(map(str, gens))} rows=(1,{q}) {' '.join(info)}", flush=True)
                    continue
                print(f"{lE2 + gain:7.2f} (E'={lE2:.2f} +{gain:.2f}) n={v} 4n={4*v} |H|={len(H)} gens={','.join(map(str, gens))} "
                      f"orbits={st.O} classes={st.C} rows=(1,{q[0]},{q[1]},{q[2]}) log10E={lE:.2f} "
                      f"log10N=({','.join(f'{x:.1f}' for x in lN)}) {' '.join(info)}", flush=True)

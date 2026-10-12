"""Screen skew Goethals-Seidel multiplier structures for composite v: exact compression test, then the model.

For each subgroup H of Z_v^* with -1 not in H and MINO..MAXO orbits, and each row-sum pattern (1; b, c, d):
1. For every prime p | v, the compressions to Z_p are listed exactly (a DP over the H-orbits, as in
   chain.py) and joined in integers. An empty join rules the structure and pattern out.
2. Survivors get the Gaussian estimate of skew_rank.py, with the low-energy gain. For composite v this
   estimate is unreliable (off by 10^2-10^5 either way), so it is only used to order the survivors; the
   level-1 counts are printed beside it.

  python3 attack/census_screen.py MINO MAXO CAP v ...     -> $HADAMARD_WORK/census_screen.jsonl
"""
import json, math, sys, time
import numpy as np
from pathlib import Path as _P; sys.path[:0] = [str(_P(__file__).resolve().parent)]
import workdir as W
import spence_lib as L
from gs_sweep_subgroups import subgroups
import chain
import skew_rank as R


def screen(v, gens, pat, maxstates=300_000):
    roles = [None] + list(pat)
    ps = [p for p in range(2, v) if v % p == 0 and all(p % r for r in range(2, int(p ** 0.5) + 1))]
    counts = {}
    for p in ps:
        c1 = chain.Chain(v, gens, [p], skew=True)
        Sd = [c1.states(r, maxstates) for r in roles]
        if any(x is None for x in Sd): continue
        U = [np.unique(K_, axis=0).astype(np.int64) for K_, _ in Sd]
        n1 = 0 if any(len(u) == 0 for u in U) else len(chain.join4([c1.paf(0, u) for u in U]))
        counts[p] = n1
        if n1 == 0: return counts, False
    return counts, True


def main():
    mino, maxo, cap = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
    out = W.WORK / 'census_screen.jsonl'
    rng = np.random.default_rng(1)
    for v in map(int, sys.argv[4:]):
        t0 = time.time(); nst = nsurv = 0
        for Hs, gens in subgroups(v).items():
            if (v - 1) in Hs: continue
            O = len({min(t * h % v for h in Hs) for t in range(v)})
            if not (mino <= O <= maxo): continue
            G = L.AbGroup([v]); st = L.Structure(G, [G.perm_from_matrix(np.array([[g]])) for g in gens])
            res = R.rank_structure(st, v, cap, rng=rng)
            for lE2, lE, q, lN in sorted(res, reverse=True):
                nst += 1
                lG = lE2 + R.GAIN * (st.C - 1)
                try:
                    counts, alive = screen(v, list(gens), q)
                except Exception as e:
                    counts, alive = {'error': repr(e)}, None
                rec = dict(v=v, gens=list(gens), H=len(Hs), O=st.O, C=st.C, pattern=list(q), log10E=round(lE, 2),
                           log10E_capped=round(lE2, 2), log10E_gen=round(lG, 2), level1=counts, alive=alive)
                with open(out, 'a') as f: f.write(json.dumps(rec) + '\n')
                if alive: nsurv += 1; print(json.dumps(rec), flush=True)
        print(f"== v={v}: {nsurv} of {nst} structure-patterns survive [{time.time()-t0:.0f}s]", flush=True)


if __name__ == '__main__':
    main()

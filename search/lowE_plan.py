"""Plan a skew_lowE.py run: Gaussian estimate of the number of solutions inside the lists for several list sizes,
from 50,000-sample statistics of low-energy candidates (skew_gen, given descent depth) for each row sum.
Roles with equal row sums share one list (as in skew_lowE.py). One class is implied by the row sums and left out.

  python3 search/lowE_plan.py v g b,c,d [--steps S]
"""
import argparse, math, sys
import numpy as np
from pathlib import Path as _P; sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import spence_lib as L
import skew_ktree as K
from skew_lowE import log10_yield


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('v', type=int); ap.add_argument('g'); ap.add_argument('pattern')
    ap.add_argument('--steps', type=int, default=4); ap.add_argument('--samples', type=int, default=50000)
    a = ap.parse_args()
    v = a.v; gens = [int(x) for x in a.g.split(',')]; pat = [int(x) for x in a.pattern.split(',')]
    G = L.AbGroup([v]); st = L.Structure(G, [G.perm_from_matrix(np.array([[g]])) for g in gens]); pm = v % 4
    rng = np.random.default_rng(5)
    rep = np.array([np.flatnonzero(st.orb == o)[0] for o in range(st.O)]); partner = st.orb[G.neg()[rep]]
    pairs = [(o, partner[o]) for o in range(1, st.O) if o < partner[o]]
    s = 1 - 2 * ((np.arange(2 ** len(pairs))[:, None] >> np.arange(len(pairs))) & 1)
    YA = np.zeros((len(s), st.O), np.int8); YA[:, 0] = 1
    for j, (o, q) in enumerate(pairs): YA[:, o] = s[:, j]; YA[:, q] = -s[:, j]
    VA, _ = K.values(st, YA, pm)
    V = {}
    for r in sorted(set(pat)):
        Y = K.gen_pool(st, r, a.samples, a.steps, rng, v); V[r], _ = K.values(st, Y, pm)
        print(f"row sum {r}: {len(Y)} distinct of {a.samples}, pool {st.pool(r):.3g}", flush=True)
    b, c, d = pat
    for nbc, nd in ((1e6, 1e6), (2e6, 8e6), (4e6, 1.6e7), (8e6, 3e7)):
        want = {}
        for r, n in ((b, nbc), (c, nbc), (d, nd)): want[r] = max(want.get(r, 0), n)
        sizes = [len(VA)] + [want[r] for r in pat]
        lg = log10_yield([VA] + [V[r] for r in pat], pm, sizes, st.csize)
        print(f"v={v} H={gens} rows=(1,{b},{c},{d}) steps={a.steps} lists B,C={nbc:.0e} D={nd:.0e} (shared by row sum): "
              f"log10 solutions in lists {lg:.2f}; distinct list volume {sum(want.values()):.2e}", flush=True)


if __name__ == '__main__':
    main()

import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)      # ../qshape.py is the base shape code
import json, sys
from qshape import kshape
from collections import defaultdict
def est(m, q, gens, cap_mult=4):
    sh = kshape([m], gens); nu = sh['nu']; ell = sh['ell']; K = sh['Ksz']; neg = sh['neg']
    cap = cap_mult * K
    # free block: dp over orbits
    dp = {0: 1}
    for o in range(nu):
        nd = defaultdict(int)
        for w, c in dp.items():
            nd[w] += 2 * c if o > 0 else c  # E: x0=+-1 (o=0 normalised x0=1 -> 1 symbol)
            for cc, mult in ((q + 1, 2), (q - 3, 2)):
                ww = w + ell[o] * cc
                if ww <= cap: nd[ww] += (mult if o > 0 else 1) * c
        dp = nd
    nfree = sum(dp.values())
    # skew block: orbit 0 fixed; pairs
    dp = {0: 1}
    for o in range(1, nu):
        if o < neg[o]:
            nd = defaultdict(int)
            for w, c in dp.items():
                nd[w] += 2 * c
                for cc, mult in ((q + 1, 2), (q - 3, 2)):
                    ww = w + 2 * ell[o] * cc
                    if ww <= cap: nd[ww] += mult * c
            dp = nd
    nskew = sum(dp.values())
    return nu, nfree, nskew
if __name__ == '__main__':
    for a in sys.argv[1:]:
        v, m, q, g = a.split(':'); print(v, m, q, g, est(int(m), int(q), json.loads(g)))

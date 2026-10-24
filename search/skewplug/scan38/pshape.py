"""PShape: SShape with H_K given by explicit permutations of the K indices (matrix actions on Z_p^k blocks, mixed K); same API as sk.SShape (hunits is None: no production)."""
import numpy as np, itertools, sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..'))      # scan38/ and the modules of search/skewplug
DATA = os.path.join(HERE, 'data')
import pk
from sk import SShape
_orig_mulperm = pk.KGrp.mulperm
def _mp(self, units):
    reg = getattr(pk.KGrp, '_preg', None)
    if reg is not None and isinstance(units, (int, np.integer)): return np.asarray(reg[int(units)])
    return _orig_mulperm(self, units)
pk.KGrp.mulperm = _mp

class PShape(SShape):
    def __init__(self, Kmods, q, perms):
        pk.KGrp._preg = [np.asarray(p, dtype=np.int64) for p in perms]
        try:
            super().__init__(Kmods, q, list(range(len(perms))))
        finally:
            pk.KGrp._preg = None
        self.hunits = None

def _tup(x):
    if isinstance(x, (list, tuple)): return tuple(_tup(y) for y in x)
    return x

def shape_from_row(row):
    """row of shapes38.json -> PShape (via enum38 group elements)"""
    import enum38 as E
    blocks = [tuple(b) for b in row['blocks']]
    mods = E.flat_mods(blocks); N = int(np.prod(mods))
    el_coords = np.array(list(itertools.product(*[range(x) for x in mods])), dtype=np.int64).reshape(N, len(mods))
    wv = np.ones(len(mods), dtype=np.int64)
    for i in range(len(mods) - 2, -1, -1): wv[i] = wv[i + 1] * mods[i + 1]
    perms = [E.perm_of(blocks, _tup(g), el_coords, wv) for g in row['gens']]
    return PShape(mods, row['q'], perms)

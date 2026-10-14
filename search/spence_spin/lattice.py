"""Joint (full-covariance) lattice gate for the full structure.  Keys: class-PAF vectors (D=27) of A (row sum 2, K-orbit signs)
and of the spin-sum of B (row sum 0, H-orbit signs).  Lambda_X = Z-span of differences of keys over the pool (sampled; a sample lattice is a
sublattice of the true one, so  'target in sample lattice'  PROVES the lattice condition; 'not in'  is only evidence until the sample lattice stabilises).
Condition: -4*1 - kA0 - kB0 in Lambda_A + Lambda_B."""
import sys, time
import numpy as np
from math import gcd
from struct573 import build

class Lat:
    """echelon basis over Z (python ints), size-reduced against pivots"""
    def __init__(self, d): self.d = d; self.rows = {}
    def _red(self, v, start):
        for c in range(start, self.d):
            r = self.rows.get(c)
            if r is not None and v[c] != 0:
                q = v[c] // r[c]
                if q: v = [vi - q * ri for ri, vi in zip(r, v)]
        return v
    def add(self, v):
        v = [int(x) for x in v]
        for c in range(self.d):
            if v[c] == 0: continue
            if c not in self.rows:
                if v[c] < 0: v = [-x for x in v]
                self.rows[c] = self._red(v, c + 1) if True else v
                # keep pivot positive
                return True
            r = self.rows[c]
            a, b = r[c], v[c]
            g, x, y = egcd(a, b)
            newp = [x * ri + y * vi for ri, vi in zip(r, v)]
            left = [(a // g) * vi - (b // g) * ri for ri, vi in zip(r, v)]
            self.rows[c] = self._red(newp, c + 1)
            v = self._red(left, c + 1)
        return False
    def contains(self, v):
        v = [int(x) for x in v]
        for c in range(self.d):
            if v[c] == 0: continue
            if c not in self.rows: return False
            r = self.rows[c]
            if v[c] % r[c]: return False
            q = v[c] // r[c]
            v = [vi - q * ri for ri, vi in zip(r, v)]
        return True
    def index_data(self):
        return [self.rows[c][c] if c in self.rows else 0 for c in range(self.d)]

def egcd(a, b):
    x0, y0, x1, y1 = 1, 0, 0, 1
    while b:
        q = a // b
        a, b = b, a - q * b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    if a < 0: a, x0, y0 = -a, -x0, -y0
    return a, x0, y0

def run(u, nsamp, seed=1):
    S = build(u)
    rng = np.random.default_rng(seed)
    XA = S.sample('A', 2, nsamp, rng); XB = S.sample('B', 0, nsamp, rng)
    KA = np.rint(S.valsA(XA)).astype(np.int64); KB = np.rint(S.valsB(XB)).astype(np.int64)
    D = S.D
    La, Lb = Lat(D), Lat(D)
    for i in range(1, nsamp):
        La.add(KA[i] - KA[0]); Lb.add(KB[i] - KB[0])
    Lab = Lat(D)
    for L in (La, Lb):
        for c, r in L.rows.items(): Lab.add(r)
    tgt = -4 - KA[0] - KB[0]
    return dict(u=u, nsamp=nsamp, rankA=len(La.rows), rankB=len(Lb.rows), rankAB=len(Lab.rows),
                diagA=La.index_data(), diagB=Lb.index_data(), diagAB=Lab.index_data(),
                target_in_A=La.contains(tgt), target_in_B=Lb.contains(tgt), target_in_AB=Lab.contains(tgt),
                detAB=int(np.prod([float(x) for x in Lab.index_data()])) if len(Lab.rows) == D else None)

if __name__ == '__main__':
    u = int(sys.argv[1]); ns = int(sys.argv[2])
    t = time.time()
    r = run(u, ns)
    print(r, '%.1fs' % (time.time() - t))

"""Exact compression test for skew Goethals-Seidel multiplier structures.

For d | v, compress each block to Z_d: x_j = sum_{t = j mod d} X(t). The compressed blocks satisfy
    sum_X |X^(chi)|^2 = 4v      for every character chi of Z_d,
and x is constant on the orbits of H acting on Z_d. For the skew block A, x_0 = 1 and x_{-j} = -x_j.
For each role we enumerate the achievable compressed vectors (sums over the H-orbits of Z_v lying over
each H-orbit of Z_d) and test whether one vector per role can meet the norm equations. A structure that
fails for some d has no solutions at all; the Gaussian model of skew_triage.py does not see this.

  python3 skew_compress.py v gens b,c,d [d ...]
"""
import itertools, math, sys
from collections import defaultdict
import numpy as np
sys.path.insert(0, '/home/user/Hadamard/attack')
import spence_lib as L


def compressed_options(st, v, d, skew=False, r=None, cap=2_000_000):
    """Set of achievable compressed vectors on Z_d (tuple over Z_d), as a dict: vector -> True."""
    orb = st.orb
    # Z_d-orbits under H (via the images of Z_v-orbits)
    # contribution of each Z_v-orbit O: a vector over Z_d (counts of O in each residue class)
    contrib = np.zeros((st.O, d), dtype=np.int64)
    for t in range(v): contrib[orb[t], t % d] += 1
    if skew:
        rep = np.array([np.flatnonzero(orb == o)[0] for o in range(st.O)])
        partner = orb[st.G.neg()[rep]]
        pairs = [(o, partner[o]) for o in range(1, st.O) if o < partner[o]]
        items = [contrib[o] - contrib[q] for o, q in pairs]          # sign +: A = +1 on o, -1 on q
        base = contrib[0].copy()
    else:
        items = [2 * contrib[o] for o in range(st.O)]                  # x = sum sign*contrib = -sum contrib + sum_{+} 2 contrib
        base = -contrib.sum(0)
    # DP over items with dedup of reachable vectors (and the row-sum constraint at the end)
    reach = {tuple(base)}
    for k, it in enumerate(items):
        if skew:
            new = {tuple(np.array(x) + it) for x in reach} | {tuple(np.array(x) - it) for x in reach}
        else:
            new = reach | {tuple(np.array(x) + it) for x in reach}
        reach = new
        if len(reach) > cap: return None
    if r is not None: reach = {x for x in reach if sum(x) == r}
    return reach


def norms(xs, d):
    """|x^(chi_k)|^2 for k = 1..d-1 (exact integers)."""
    X = np.array(sorted(xs), dtype=np.float64)
    F = np.fft.fft(X, axis=1)
    N = np.rint((F.real ** 2 + F.imag ** 2)[:, 1:]).astype(np.int64)
    return {tuple(row) for row in N}


def feasible(v, gens, pattern, d, verbose=True):
    G = L.AbGroup([v]); st = L.Structure(G, [G.perm_from_matrix(np.array([[g]])) for g in gens])
    opts = [compressed_options(st, v, d, skew=True)] + [compressed_options(st, v, d, r=r) for r in pattern]
    if any(o is None for o in opts): return None
    Ns = [norms(o, d) for o in opts]
    target = 4 * v
    left = defaultdict(int)
    for a in Ns[0]:
        for b in Ns[1]: left[tuple(x + y for x, y in zip(a, b))] += 1
    hits = 0
    for c in Ns[2]:
        for e in Ns[3]:
            need = tuple(target - x - y for x, y in zip(c, e))
            if need in left: hits += 1
    if verbose: print(f"v={v} gens={gens} rows=(1,{pattern}) d={d}: norm vectors {[len(n) for n in Ns]}, feasible combinations {hits}")
    return hits > 0


if __name__ == '__main__':
    v = int(sys.argv[1]); gens = [int(x) for x in sys.argv[2].split(',')]; pat = [int(x) for x in sys.argv[3].split(',')]
    ds = [int(x) for x in sys.argv[4:]] or [p for p in range(2, v) if v % p == 0 and all(p % q for q in range(2, int(p ** 0.5) + 1))]
    for d in ds: feasible(v, gens, pat, d)

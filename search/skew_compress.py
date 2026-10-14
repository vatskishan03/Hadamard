"""Exact compression test for skew Goethals-Seidel multiplier structures.

For d | v, compress each block to Z_d: x_j = sum_{t = j mod d} X(t). The compressed blocks are complementary:
    sum_X PAF_x(k) = 4v [k = 0]      for k in Z_d   (integers; equivalently sum_X |X^(chi)|^2 = 4v for all chi),
and x is constant on the orbits of H acting on Z_d. For the skew block A, x_0 = 1 and x_{-j} = -x_j.
For each role we enumerate the achievable compressed vectors (sums over the H-orbits of Z_v lying over
each H-orbit of Z_d) and test whether one vector per role can meet the norm equations. A structure that
fails for some d has no solutions at all; the Gaussian model of skew_triage.py does not see this.

  python3 skew_compress.py v gens b,c,d [d ...]
"""
import itertools, math, sys
from collections import defaultdict
import numpy as np
from pathlib import Path as _P; sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
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


def pafs(xs, d):
    """Exact integer periodic autocorrelations (shifts 0..d//2) of compressed vectors on Z_d.
    Summed over the four blocks they must equal (4v, 0, ..., 0): the compression of a complementary family
    is complementary. This is equivalent to the norm equations at every character, but it stays in the integers
    (|x^(chi)|^2 is not a rational integer when d >= 5, so rounded norms are not a valid test)."""
    X = np.array(sorted(xs), dtype=np.int64)
    P = np.stack([(X * np.roll(X, -k, axis=1)).sum(1) for k in range(d // 2 + 1)], axis=1)
    return {tuple(row) for row in P}


def target(v, d):
    return tuple([4 * v] + [0] * (d // 2))


def feasible(v, gens, pattern, d, verbose=True):
    G = L.AbGroup([v]); st = L.Structure(G, [G.perm_from_matrix(np.array([[g]])) for g in gens])
    opts = [compressed_options(st, v, d, skew=True)] + [compressed_options(st, v, d, r=r) for r in pattern]
    if any(o is None for o in opts): return None
    Ns = [pafs(o, d) for o in opts]
    T = target(v, d)
    left = defaultdict(int)
    for a in Ns[0]:
        for b in Ns[1]: left[tuple(x + y for x, y in zip(a, b))] += 1
    hits = 0
    for c in Ns[2]:
        for e in Ns[3]:
            need = tuple(t - x - y for t, x, y in zip(T, c, e))
            if need in left: hits += left[need]
    if verbose: print(f"v={v} gens={gens} rows=(1,{pattern}) d={d}: compressed autocorrelation vectors {[len(n) for n in Ns]}, feasible combinations {hits}")
    return hits > 0


def role_filters(st, v, d, pattern, cap=2_000_000):
    """For each role (A, then the row sums in pattern), the set of compressed norm vectors that occur in at
    least one feasible combination; None if the enumeration is too large."""
    opts = [compressed_options(st, v, d, skew=True, cap=cap)] + [compressed_options(st, v, d, r=r, cap=cap) for r in pattern]
    if any(o is None for o in opts): return None
    Ns = [sorted(pafs(o, d)) for o in opts]
    T = target(v, d)
    left = defaultdict(list)
    for a in Ns[0]:
        for b in Ns[1]: left[tuple(x + y for x, y in zip(a, b))].append((a, b))
    keep = [set(), set(), set(), set()]
    for c in Ns[2]:
        for e in Ns[3]:
            need = tuple(t - x - y for t, x, y in zip(T, c, e))
            for a, b in left.get(need, ()):
                keep[0].add(a); keep[1].add(b); keep[2].add(c); keep[3].add(e)
    return keep


def norm_rows(Y, st, d):
    """Exact compressed autocorrelation vectors (shifts 0..d//2) of orbit sign vectors Y (rows); same keys as pafs()."""
    contrib = np.zeros((st.O, d), dtype=np.int64)
    for t in range(st.G.v): contrib[st.orb[t], t % d] += 1
    X = np.rint(Y).astype(np.int64) @ contrib
    return np.stack([(X * np.roll(X, -k, axis=1)).sum(1) for k in range(d // 2 + 1)], axis=1)

if __name__ == '__main__':
    v = int(sys.argv[1]); gens = [int(x) for x in sys.argv[2].split(',')]; pat = [int(x) for x in sys.argv[3].split(',')]
    ds = [int(x) for x in sys.argv[4:]] or [p for p in range(2, v) if v % p == 0 and all(p % q for q in range(2, int(p ** 0.5) + 1))]
    for d in ds: feasible(v, gens, pat, d)

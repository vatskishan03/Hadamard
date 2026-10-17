"""Paley plus three blocks: skew Goethals-Seidel quadruples of prime length v = 3 (mod 4) whose first block is Paley's
(-1 on the nonzero squares, PAF -1 off zero) and whose other three blocks are unions of orbits of a subgroup H of the
squares, of odd order d, with PAF sum +1 off zero (the difference families of Djokovic-Kotsireas, arXiv:1511.08734).

    python3 search/skew_paley3.py exhaust v d [limit]    exhaustive search: PAF as quadratic forms in the orbit signs,
                                                         PSD bound, first block up to the action of the squares, hashed
                                                         join of the third block
    python3 search/skew_paley3.py sweep lo hi             every census-undecided prime v = 3 (mod 4) and every d with
                                                         lo <= (v - 1)/d <= hi
    python3 search/skew_paley3.py estimate v/d ...        Gaussian count of solutions (ordered triples, sign of the 0
                                                         orbit fixed), the row-sum direction removed

search/NOTES.md, section 48.
"""
import itertools, math, sys, time
from math import comb, log10
from pathlib import Path
import numpy as np


def run(v, d, limit=1, log=print):
    g = next(x for x in range(2, v) if len({pow(x, k, v) for k in range(v - 1)}) == v - 1)
    h = pow(g, (v - 1) // d, v)
    H = sorted({pow(h, k, v) for k in range(d)})
    m = (v - 1) // d
    orb_of = np.zeros(v, dtype=np.int64)
    orbits = [[0]]
    reps = []
    for j in range(m):                       # orbit j+1 = g^j H
        r = pow(g, j, v)
        o = sorted({r * x % v for x in H})
        orbits.append(o); reps.append(r)
        orb_of[o] = j + 1
    nb = m + 1
    # N[k][i][j] = #{x in O_i : x + s_k in O_j}, s_k = reps[k]
    N = np.zeros((m, nb, nb), dtype=np.int64)
    for k, s in enumerate(reps):
        for i, o in enumerate(orbits):
            for x in o:
                N[k, i, orb_of[(x + s) % v]] += 1
    # Gauss periods for PSD: G[i][k] = sum_{x in O_i} w^(s_k x)
    w = np.exp(2j * np.pi / v)
    G = np.array([[sum(w ** (s * x % v) for x in o) for s in reps] for o in orbits])
    # all sign vectors with c_0 = +1
    M = 1 << m
    C = np.ones((M, nb), dtype=np.int8)
    idx = np.arange(M)
    for b in range(m):
        C[(idx >> b) & 1 == 1, b + 1] = -1
    Cf = C.astype(np.float64)
    psd = np.abs(Cf @ G) ** 2
    ok = np.all(psd <= 3 * v - 1 + 1e-6, axis=1)
    sizes = np.array([len(o) for o in orbits])
    rs = (C.astype(np.int64) @ sizes)
    C, Cf, psd, rs = C[ok], Cf[ok], psd[ok], rs[ok]
    log(f"v {v} |H| {d} orbits {m}: {M} sign vectors, {len(C)} pass the PSD bound")
    # PAF on the m shift representatives
    P = np.empty((len(C), m), dtype=np.int64)
    Ci = C.astype(np.int64)
    for k in range(m):
        P[:, k] = np.einsum("ni,ij,nj->n", Ci, N[k], Ci, optimize=True)
    # first-block canonical forms under multiplication by g (g^j H = orbit j+1): rotation of orbits 1..m by an even
    # amount keeps Paley invariant (squares); g^2 rotates the orbit indices by 2 (mod m)
    def canon_mask():
        keys = C[:, 1:].astype(np.int8)
        best = np.ones(len(C), dtype=bool)
        pw = (1 << np.arange(m, dtype=np.int64))
        code = ((keys < 0).astype(np.int64) @ pw)
        for t in range(2, m, 2):
            rot = np.roll(keys, t, axis=1)
            c2 = ((rot < 0).astype(np.int64) @ pw)
            best &= code <= c2
        return best
    canon = canon_mask()
    absr = np.abs(rs)
    target = 4 * v - 1
    vals = sorted(set(absr.tolist()))
    rng = np.random.default_rng(1)
    HV = rng.integers(1, 2 ** 31, size=m)
    found = []
    for r1 in vals:
        for r2 in vals:
            if r2 < r1: continue
            r3sq = target - r1 * r1 - r2 * r2
            if r3sq < r2 * r2: continue
            r3 = int(round(r3sq ** 0.5))
            if r3 * r3 != r3sq or r3 not in vals: continue
            L1 = np.flatnonzero((absr == r1) & canon)
            L2 = np.flatnonzero(absr == r2)
            L3 = np.flatnonzero(absr == r3)
            keys3 = P[L3] @ HV
            order = np.argsort(keys3); ks = keys3[order]
            t0 = time.time()
            for a in L1:
                Tt = 1 - P[a][None, :] - P[L2]
                kk = Tt @ HV
                pos = np.searchsorted(ks, kk)
                pos[pos >= len(ks)] = 0
                hit = np.flatnonzero(ks[pos] == kk)
                for hi in hit:
                    b = L2[hi]
                    c = L3[order[pos[hi]]]
                    if np.array_equal(P[a] + P[b] + P[c], np.ones(m, dtype=np.int64)):
                        sets = [sorted(x for i in range(nb) if C[q][i] < 0 for x in orbits[i]) for q in (a, b, c)]
                        found.append(sets)
                        log(f"HIT v {v} |H| {d} rows {r1} {r2} {r3}: {sets}")
                        if len(found) >= limit:
                            return H, found
            log(f"  rows {r1} {r2} {r3}: |L1| {len(L1)} |L2| {len(L2)} |L3| {len(L3)} "
                f"joined in {time.time() - t0:.1f}s")
    return H, found


def est(v, d, nsamp=1500, seed=0):
    g = next(x for x in range(2, v) if len({pow(x, k, v) for k in range(v - 1)}) == v - 1)
    h = pow(g, (v - 1) // d, v); H = sorted({pow(h, k, v) for k in range(d)})
    m = (v - 1) // d
    orbits = [[0]] + [sorted({pow(g, j, v) * x % v for x in H}) for j in range(m)]
    reps = [pow(g, j, v) for j in range(m // 2)]
    rng = np.random.default_rng(seed)
    def sample(j, n):
        X = np.ones((n, v))
        for t in range(n):
            for o in rng.choice(m, j, replace=False): X[t, orbits[o + 1]] = -1
        return np.rint(np.fft.ifft(np.abs(np.fft.fft(X, axis=1)) ** 2, axis=1).real)[:, reps]
    tgt = 4 * v - 1
    rowj = {}
    for j in range(m + 1):
        r = abs(1 + d * (m - 2 * j)); rowj.setdefault(r, []).append(j)
    vals = sorted(rowj); trip = []
    for a in vals:
        for b in vals:
            if b < a: continue
            c2 = tgt - a*a - b*b
            if c2 < b*b: continue
            c = int(round(c2 ** .5))
            if c*c == c2 and c in rowj: trip.append((a, b, c))
    k = m // 2
    stats = {}
    tot = 0.0
    for a, b, c in trip:
        for r in (a, b, c):
            if r not in stats:
                Ps = np.vstack([sample(j, max(50, nsamp // len(rowj[r]))) for j in rowj[r]])
                stats[r] = (Ps.mean(0), np.cov(Ps.T).reshape(k, k), sum(comb(m, j) for j in rowj[r]))
        mu = stats[a][0] + stats[b][0] + stats[c][0]
        Sg = stats[a][1] + stats[b][1] + stats[c][1]
        x = np.ones(k) - mu
        if k > 1:
            Sr, xr = Sg[:-1, :-1], x[:-1]
            sign, logdet = np.linalg.slogdet(Sr)
            lp = (-0.5 * xr @ np.linalg.solve(Sr, xr) - 0.5 * logdet - (k - 1) / 2 * np.log(2 * np.pi)
                  + (k - 1) * np.log(4))      # PAF values lie on a lattice of spacing 4
        else:
            lp = 0.0
        perms = 6 if len({a, b, c}) == 3 else (3 if len({a, b, c}) == 2 else 1)
        cnt = sum(math.log(stats[r][2]) for r in (a, b, c)) + math.log(perms)
        tot += math.exp(min(cnt + lp, 700))
    return trip, (log10(tot) if tot > 0 else float('-inf'))


def undecided():
    p = Path(__file__).resolve().parent / "data" / "skew_census_undecided.txt"
    return [int(x) for x in p.read_text().splitlines() if x.strip() and not x.startswith("#")]


def is_prime(n):
    return n > 1 and all(n % q for q in range(2, math.isqrt(n) + 1))


if __name__ == "__main__":
    cmd = sys.argv[1]
    log = lambda s: print(s, flush=True)
    if cmd == "exhaust":
        v, d = int(sys.argv[2]), int(sys.argv[3])
        H, f = run(v, d, int(sys.argv[4]) if len(sys.argv) > 4 else 1, log=log)
        print("solutions found:", len(f))
    elif cmd == "sweep":
        lo, hi = int(sys.argv[2]), int(sys.argv[3])
        jobs = sorted(((v - 1) // d, v, d) for v in undecided() if is_prime(v) and v % 4 == 3
                      for d in range(3, v, 2) if (v - 1) % (2 * d) == 0 and lo <= (v - 1) // d <= hi)
        for orb, v, d in jobs:
            t0 = time.time()
            H, f = run(v, d, 1, log=log)
            print(f"RESULT v {v} |H| {d} orbits {orb}: {'HIT' if f else 'none'} ({time.time() - t0:.0f}s)", flush=True)
    elif cmd == "estimate":
        for vd in sys.argv[2:]:
            v, d = map(int, vd.split("/"))
            trip, le = est(v, d)
            print(f"v {v} |H| {d} orbits {(v - 1) // d}: row-sum triples {trip}, log10 E = {le:.2f}", flush=True)

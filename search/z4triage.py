"""Triage of the bordered families over Q = Z_4 (NOTES section 22) for the FRONTIER orders n = 3 (mod 4).

G = Z_N (N = n - 3, cyclic), K = 4 Z_N, T = {1, chi, chi^3} ("013", Gamma = 12, 4, -4, 4) or T = {chi, chi^2, chi^3}
("123", Gamma = 12, -4, -4, -4); targets sum_i PAF_i(g) = -Gamma(g mod 4) for g != 0, border of size 12
(search/border_probe.py). The sequences are invariant under a multiplier group H <= Z_N^* (every multiplier maps
g mod 4 to +-g mod 4 and Gamma(b) = Gamma(-b), so the targets are H-invariant).

For each H with between OMIN and OMAX orbits:
1. stage 1 on Q = Z_4, exactly: the compressions y_i(b) (sums over the cosets of K) of H-invariant sequences with row
   sums r_i must satisfy sum_i |y^_i(j)|^2 = 12 for j in T and 4n for j not in T. Each role's achievable compressions
   are listed by a dynamic program over the orbits, and the four are joined.
2. If stage 1 is not empty: the Gaussian estimate of the number of solutions (full covariance of the class values of
   sampled members, lattice spacing from the sampled values), per row-sum pattern.

  python3 search/z4triage.py n1,n2,.. [--omin 12] [--omax 40] [--samples 20000]   -> $HADAMARD_WORK/z4triage.jsonl
"""
import argparse, itertools, json, math, sys, time
from math import gcd, log10
import numpy as np
from pathlib import Path as _P
sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
import spence_lib as L
from gs_sweep_subgroups import subgroups

GAMMA = {'013': [12, 4, -4, 4], '123': [12, -4, -4, -4]}
TSET = {'013': (0, 1, 3), '123': (1, 2, 3)}


def row_patterns(total):
    """Unordered (r_0 <= .. <= r_3), r_i >= 0 even, sum r_i^2 = total (N even: r_i even)."""
    out = []
    lim = math.isqrt(total)
    for a in range(0, lim + 1, 2):
        for b in range(a, lim + 1, 2):
            for c in range(b, lim + 1, 2):
                d2 = total - a * a - b * b - c * c
                if d2 < c * c: continue
                d = math.isqrt(d2)
                if d * d == d2 and d % 2 == 0: out.append((a, b, c, d))
    return out


def coset_sums(st, which):
    """Achievable sums over the orbits contained in the cosets `which` (a tuple of residues mod 4) of H-invariant sign
    vectors: dict value -> count (values as tuples, one entry per residue in `which`)."""
    states = {tuple(0 for _ in which): 1}
    for o in range(st.O):
        el = np.flatnonzero(st.orb == o)
        c = [int((el % 4 == b).sum()) for b in which]
        if not any(c): continue
        new = {}
        for y, cnt in states.items():
            for sg in (1, -1):
                y2 = tuple(y[i] + sg * c[i] for i in range(len(which)))
                new[y2] = new.get(y2, 0) + cnt
        states = new
    return states


def spectra(st, N, r):
    """Achievable spectra (|y^(0)|^2, |y^(1)|^2, |y^(2)|^2) of the Z_4-compressions of H-invariant sign vectors with
    row sum +-r (|y^(3)| = |y^(1)|)."""
    if not hasattr(st, '_cs'):
        st._cs = (coset_sums(st, (0,)), coset_sums(st, (2,)), coset_sums(st, (1, 3)))
    E0, E2, E13 = st._cs
    out = set()
    for (y0,) in E0:
        for (y2,) in E2:
            for (y1, y3) in E13:
                s0 = y0 + y1 + y2 + y3
                if abs(s0) != r: continue
                out.add((s0 * s0, (y0 - y2) ** 2 + (y1 - y3) ** 2, (y0 - y1 + y2 - y3) ** 2))
    return out


def stage1(st, N, n, T, rows):
    """Number of quadruples of achievable compression spectra (one per role) with sum_i |y^_i(j)|^2 = 12 for j in T
    and 4n for j not in T (capped at 10^6); 0 means stage 1 is empty."""
    target = tuple(12 if j in T else 4 * n for j in range(3))      # j = 0, 1, 2 (j = 3 equals j = 1)
    L_ = [sorted(spectra(st, N, r)) for r in rows]
    AB = {}
    for a in L_[0]:
        for b in L_[1]:
            sab = tuple(x + y for x, y in zip(a, b))
            if all(sab[j] <= target[j] for j in range(3)): AB[sab] = AB.get(sab, 0) + 1
    Dset = set(L_[3]); found = 0
    for sab, cnt in AB.items():
        for c in L_[2]:
            need = tuple(target[j] - sab[j] - c[j] for j in range(3))
            if need in Dset:
                found += cnt
                if found > 10 ** 6: return found
    return found


def gauss(st, rows, gam, samples, rng):
    """log10 of the expected number of solutions for one row pattern (full covariance; lattice from the samples)."""
    keys, logs = [], 0.0
    for r in rows:
        N_ = st.pool(r)
        if N_ <= 0: return -math.inf
        logs += log10(N_)
        keys.append(st.class_values(st.sample(r, samples, rng)))
    tgt = np.array([-gam[int(st.reps[c]) % 4] for c in range(st.C)], dtype=float)
    mu = sum(k.mean(0) for k in keys); Sig = sum(np.cov(k.T) for k in keys)
    w_, U = np.linalg.eigh(Sig); keep = w_ > 1e-6 * w_.max()
    off = U[:, ~keep].T @ (tgt - mu)
    if keep.sum() < len(w_) and np.abs(off).max() > 1e-6: return -math.inf
    z = U[:, keep].T @ (tgt - mu)
    ld = -0.5 * (keep.sum() * math.log(2 * math.pi) + np.log(w_[keep]).sum()) - 0.5 * float((z * z / w_[keep]).sum())
    lat = 0.0
    if keep.all():
        for c in range(st.C):
            g = 0
            for k in keys:
                u = np.unique(np.rint(k[:, c]).astype(np.int64))
                if len(u) > 1: g = gcd(g, int(np.gcd.reduce(np.diff(u))))
            if g > 1: lat += log10(g)
    return logs + ld / math.log(10) + lat


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('ns'); ap.add_argument('--omin', type=int, default=12); ap.add_argument('--omax', type=int, default=40)
    ap.add_argument('--samples', type=int, default=20000); ap.add_argument('--out', default='z4triage.jsonl')
    a = ap.parse_args()
    out = W.WORK / a.out
    rng = np.random.default_rng(1)
    for n in [int(x) for x in a.ns.split(',')]:
        N = n - 3
        G = L.AbGroup([N])
        P = lambda g: G.perm_from_matrix(np.array([[g]]))
        for Hs, gens in sorted(subgroups(N).items(), key=lambda t: -len(t[0])):
            if (N - 1) // len(Hs) + 1 > 3 * a.omax: continue
            t0 = time.time()
            st = L.Structure(G, [P(g) for g in gens])
            if not (a.omin <= st.O <= a.omax): continue
            for tk in ('013', '123'):
                T = TSET[tk]; gam = GAMMA[tk]
                total = 12 if tk == '013' else 4 * n
                for rows in row_patterns(total):
                    rec = dict(n=n, N=N, H=len(Hs), gens=list(gens), O=int(st.O), C=int(st.C), T=tk, rows=list(rows))
                    s1 = stage1(st, N, n, T, rows)
                    rec['stage1'] = s1
                    if s1:
                        rec['logE'] = round(gauss(st, rows, gam, a.samples, rng), 2)
                        rec['pools'] = [round(log10(max(st.pool(r), 1)), 1) for r in rows]
                    rec['sec'] = round(time.time() - t0, 1)
                    print(json.dumps(rec), flush=True)
                    with open(out, 'a') as f: f.write(json.dumps(rec) + '\n')


if __name__ == '__main__':
    main()

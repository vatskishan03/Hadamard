r"""Gaussian estimate and tabu search for the tri-skew structure (A, A o mu, A o mu^2, D) of tri_skew.py.

A is skew and H-invariant (one sign per pair of opposite H-orbits), D is K-invariant (K = <H, mu>) with row sum
d = sqrt(4v - 3), and the condition is T_A + P_D = 0 on the <K, -1>-classes, T_A(s) = sum_k P_A(mu^k s).  The tabu
search is spin_ls2.c with the roles exchanged: its "B" block is D (row-sum preserving swaps of equal K-orbits) and
its "A" block is A over the H-orbits, with the quadratic form of T_A.

  python3 search/tri_ls.py rank v                    every (H, mu) of tri_skew.py --sweep with at most 400 H-orbits
  python3 search/tri_ls.py search v hgens mu [--restarts R] [--steps S] [--seed X]   (hgens comma-separated)
"""
import argparse
import json
import math
import os
import subprocess
import sys
import time
from math import gcd, log10, pi

import numpy as np

sys.path[:0] = [os.path.dirname(os.path.abspath(__file__)),
                os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'verify')]
import workdir as W  # noqa: E402
import spence_lib as L  # noqa: E402
import verify  # noqa: E402


def closure(gens, v):
    S = {1}
    while True:
        T = S | {s * g % v for s in S for g in gens}
        if T == S:
            return frozenset(S)
        S = T


class Tri:
    def __init__(self, v, hgens, mu):
        self.v, self.mu = v, mu
        G = L.AbGroup([v])
        self.G = G
        P = lambda g: G.perm_from_matrix(np.array([[g]]))
        self.stH = L.Structure(G, [P(g) for g in hgens])
        self.stK = L.Structure(G, [P(g) for g in hgens] + [P(mu)])
        self.H, self.K = closure(list(hgens), v), closure(list(hgens) + [mu], v)
        assert (v - 1) not in self.K, '-1 in K'
        assert len(self.K) == 3 * len(self.H), 'mu is not of order 3 modulo H'
        stH, stK = self.stH, self.stK
        rep = np.array([np.flatnonzero(stH.orb == o)[0] for o in range(stH.O)])
        partner = stH.orb[G.neg()[rep]]
        self.pairs = [(o, int(partner[o])) for o in range(1, stH.O) if o < partner[o]]
        self.SA = np.zeros((stK.C, stH.O, stH.O))
        for c, s in enumerate(stK.reps):
            for k in range(3):
                self.SA[c] += stH.S[stH.cls[s * pow(mu, k, v) % v]]
        self.SA = (self.SA + self.SA.transpose(0, 2, 1)) / 2
        self.MD = (stK.S + stK.S.transpose(0, 2, 1)) / 2
        self.d = math.isqrt(4 * v - 3)
        assert self.d ** 2 == 4 * v - 3, 'v is not k^2 + k + 1'

    def A_vectors(self, s):
        Y = np.zeros((len(s), self.stH.O))
        Y[:, 0] = 1
        for j, (o, q) in enumerate(self.pairs):
            Y[:, o] = s[:, j]
            Y[:, q] = -s[:, j]
        return Y

    def estimate(self, samples=3000, rng=None):
        rng = rng or np.random.default_rng(1)
        nA = len(self.pairs)
        s = np.where(rng.random((samples, nA)) < 0.5, 1.0, -1.0)
        YA = self.A_vectors(s)
        VA = np.einsum('bo,cop,bp->bc', YA, self.SA, YA, optimize=True)
        ND = self.stK.pool(self.d)
        if ND == 0:
            return None
        YD = self.stK.sample(self.d, samples, rng)
        VD = np.einsum('bo,cop,bp->bc', YD, self.MD, YD, optimize=True)
        V = VA + VD
        lp = []
        for c in range(self.stK.C):
            var, m = V[:, c].var(), V[:, c].mean()
            if var < 1e-9:
                if abs(m) > 1e-6:
                    return None
                continue
            lp.append(-0.5 * log10(2 * pi * var) - m ** 2 / (2 * var) / math.log(10) + log10(4))
        # the class sums are tied by one linear relation (the row sums), so drop the least informative class
        lp.sort()
        return nA * log10(2) + log10(ND) + sum(lp[1:]), nA * log10(2), log10(ND)

    def sequences(self, A, D):
        a = A.astype(np.int64)[self.stH.orb]
        d = D.astype(np.int64)[self.stK.orb]
        t = np.arange(self.v)
        return [a, a[t * self.mu % self.v], a[t * self.mu * self.mu % self.v], d]


def odd_subgroups(v):
    units = [u for u in range(1, v) if gcd(u, v) == 1]
    def order(u):
        x, e = u, 1
        while x != 1:
            x, e = x * u % v, e + 1
        return e
    odd = [u for u in units if order(u) % 2 == 1]
    cyc = {}
    for u in odd:
        cyc.setdefault(closure([u], v), u)
    subs = {c: [g] for c, g in cyc.items()}
    grow = True
    while grow:
        grow = False
        for c1, g1 in list(subs.items()):
            for c2, g2 in cyc.items():
                if c2 <= c1:
                    continue
                c = closure(g1 + [g2], v)
                if c not in subs:
                    subs[c] = g1 + [g2]
                    grow = True
    return subs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=('rank', 'search'))
    ap.add_argument('v', type=int)
    ap.add_argument('h', nargs='?', default='1')
    ap.add_argument('mu', nargs='?', type=int)
    ap.add_argument('--restarts', type=int, default=10 ** 9)
    ap.add_argument('--steps', type=int, default=3000)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--maxorbits', type=int, default=400)
    a = ap.parse_args()
    v = a.v
    if a.cmd == 'rank':
        rng = np.random.default_rng(1)
        for Hs, gens in sorted(odd_subgroups(v).items(), key=lambda t: -len(t[0])):
            if (v - 1) // len(Hs) + 1 > a.maxorbits:
                continue
            seen = set()
            for mu in range(2, v):
                if gcd(mu, v) != 1 or mu in Hs or pow(mu, 3, v) not in Hs:
                    continue
                Ks = closure(list(gens) + [mu], v)
                if (v - 1) in Ks or Ks in seen:
                    continue
                seen.add(Ks)
                try:
                    t = Tri(v, list(gens), mu)
                except AssertionError:
                    continue
                e = t.estimate(rng=rng)
                if e is None:
                    print(f'v={v} |H|={len(Hs)} gens={gens} mu={mu} OH={t.stH.O} OK={t.stK.O} classes={t.stK.C}: no D', flush=True)
                    continue
                print(f'v={v} |H|={len(Hs)} gens={gens} mu={mu} OH={t.stH.O} OK={t.stK.O} classes={t.stK.C} '
                      f'pairs={len(t.pairs)} log10E={e[0]:.2f} log10NA={e[1]:.1f} log10ND={e[2]:.1f}', flush=True)
        return
    t = Tri(v, [int(x) for x in a.h.split(',')], a.mu)
    stH, stK = t.stH, t.stK
    print(f'v={v} |H|={len(t.H)} |K|={len(t.K)} OH={stH.O} OK={stK.O} classesK={stK.C} pairs={len(t.pairs)} d={t.d}',
          flush=True)
    MD2 = np.rint(2 * t.MD).astype(np.int32)
    SA2 = np.rint(2 * t.SA).astype(np.int32)
    assert np.allclose(MD2, 2 * t.MD) and np.allclose(SA2, 2 * t.SA)
    D0s = stK.sample(t.d, 4096, np.random.default_rng(a.seed)).astype(np.int32)
    path = os.path.join(W.tmpdir(), f'tri_{v}_{a.seed}.bin')
    with open(path, 'wb') as f:      # spin_ls2 input with the roles exchanged: its B is D, its A is A
        np.array([stK.O, stK.C, stH.O, len(t.pairs)], np.int32).tofile(f)
        stK.csize.astype(np.int32).tofile(f)
        stK.sizes.astype(np.int32).tofile(f)
        MD2.tofile(f)
        SA2.tofile(f)
        np.array(t.pairs, np.int32).ravel().tofile(f)
        D0s[0].tofile(f)
        np.array([len(D0s)], np.int32).tofile(f)
        D0s.tofile(f)
    t0 = time.time()
    out = subprocess.run([W.bin(os.environ.get('SPIN_LS', 'spin_ls2')), path, str(a.restarts), str(a.steps), str(a.seed)],
                         stdout=subprocess.PIPE, text=True)          # progress lines go straight to stderr
    os.remove(path)
    print(f'[{time.time() - t0:.0f}s]', flush=True)
    sol = [ln for ln in out.stdout.splitlines() if ln.startswith('SOL')]
    if not sol:
        return
    tok = sol[0].split()
    ia, ib = tok.index('A'), tok.index('B')
    A = np.array(tok[ia + 1:ib], dtype=float)
    D = np.array(tok[ib + 1:], dtype=float)
    seqs = t.sequences(A, D)
    tot = verify.paf_sum(seqs)
    ok = tot[0] == 4 * v and not tot[1:].any()
    Hm = verify.goethals_seidel(*map(verify.circulant, seqs)).astype(np.float64)
    had = np.array_equal(Hm @ Hm.T, 4 * v * np.eye(4 * v))
    skew = np.array_equal(Hm + Hm.T, 2 * np.eye(4 * v))
    print(f'FOUND: complementary={ok} Hadamard={had} skew={skew}', flush=True)
    if ok and had and skew:
        with open(W.hits(), 'a') as f:
            f.write(json.dumps(dict(n=v, gens=sorted(int(x) for x in a.h.split(',')), tri_multiplier=a.mu,
                                    pattern=[int(x.sum()) for x in seqs], hadamard=True, skew=True,
                                    method=f'tri_ls seed={a.seed} steps={a.steps}',
                                    seqs=[''.join('+' if x > 0 else '-' for x in s) for s in seqs])) + '\n')


if __name__ == '__main__':
    main()

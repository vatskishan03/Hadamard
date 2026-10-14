"""Twin skew structure (d = 2): blocks (A, A o mu, C, D) for a skew Goethals-Seidel family over Z_v.

A is skew-type and invariant under H (-1 not in H); the second block is A(mu t), also of skew type; C and D
are invariant under K = <H, mu>. Then P_{A o mu}(s) = P_A(mu s), and with T_A(s) = P_A(s) + P_A(mu s)
the family is complementary iff
    T_A(s) + P_C(s) + P_D(s) = 0     on the orbits of <H, mu, -1> on the nonzero shifts.
mu only needs mu^2 in <H, -1>; it may generate -1, because A alone must avoid it. The row sums are
(1, 1, c, d) with c^2 + d^2 = 4v - 2.

  python3 twin_skew.py rank v ...
  python3 twin_skew.py run v hgens mu c,d [cap]
"""
import json, math, os, struct, subprocess, sys, tempfile, time
from math import gcd, log10, pi
import numpy as np
from pathlib import Path as _P; sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
import spence_lib as L
from gs_sweep_subgroups import subgroups
from skew_ktree import pool_exact
import verify
HERE = str(W.HERE)


def closure(gens, v):
    S = {1}
    while True:
        T = S | {s * g % v for s in S for g in gens}
        if T == S: return frozenset(S)
        S = T


class Twin:
    def __init__(self, v, hgens, mu):
        self.v, self.mu = v, mu
        G = L.AbGroup([v]); self.G = G
        P = lambda g: G.perm_from_matrix(np.array([[g]]))
        self.stH = L.Structure(G, [P(g) for g in hgens])
        self.stK = L.Structure(G, [P(g) for g in hgens] + [P(mu)])
        stH, stK = self.stH, self.stK
        rep = np.array([np.flatnonzero(stH.orb == o)[0] for o in range(stH.O)])
        partner = stH.orb[G.neg()[rep]]
        assert not np.any(partner[1:] == np.arange(1, stH.O))
        self.pairs = [(o, partner[o]) for o in range(1, stH.O) if o < partner[o]]
        # T_A on K-class c: P_A at the H-classes of s_c and mu s_c
        self.MA = np.zeros((stK.C, stH.O, stH.O))
        for c, s in enumerate(stK.reps):
            self.MA[c] += stH.S[stH.cls[s]] + stH.S[stH.cls[s * mu % v]]
        self.MA = (self.MA + self.MA.transpose(0, 2, 1)) / 2

    def A_vectors(self, sg):
        Y = np.zeros((len(sg), self.stH.O)); Y[:, 0] = 1
        for j, (o, q) in enumerate(self.pairs): Y[:, o] = sg[:, j]; Y[:, q] = -sg[:, j]
        return Y

    def TA(self, YA):
        return np.einsum("bo,cop,bp->bc", YA, self.MA, YA, optimize=True)

    def estimate(self, c, d, samples=3000, rng=None):
        stK = self.stK
        NC, ND = stK.pool(c), stK.pool(d)
        if NC == 0 or ND == 0: return None
        sg = np.where(rng.random((samples, len(self.pairs))) < 0.5, 1.0, -1.0)
        V = self.TA(self.A_vectors(sg)) + stK.class_values(stK.sample(c, samples, rng)) + stK.class_values(stK.sample(d, samples, rng))
        lp = []
        for k in range(stK.C):
            var, m = V[:, k].var(), V[:, k].mean()
            if var < 1e-9:
                if abs(m) > 1e-6: return None
                continue
            lp.append((stK.csize[k] * math.sqrt(var), -0.5 * log10(2 * pi * var) - m * m / (2 * var) / math.log(10) + log10(4)))
        lp.sort()
        lNA = len(self.pairs) * log10(2)
        return lNA + log10(NC) + log10(ND) + sum(t for _, t in lp[:-1]), lNA, log10(NC), log10(ND)


def two_squares(N):
    out = []
    for c in range(1, int(math.isqrt(N)) + 1, 2):
        d2 = N - c * c
        if d2 < c * c: break
        d = math.isqrt(d2)
        if d * d == d2 and d % 2: out.append((c, d))
    return out


def mus(v, H):
    L_ = frozenset(H) | frozenset((v - h) % v for h in H)
    out, seen = [], set()
    for u in range(2, v):
        if gcd(u, v) != 1 or u in L_ or (u * u) % v not in L_: continue
        key = frozenset(u * x % v for x in L_)
        if key in seen: continue
        seen.add(key); out.append(u)
    return out


if __name__ == '__main__':
    if sys.argv[1] == 'rank':
        rng = np.random.default_rng(1)
        for v in map(int, sys.argv[2:]):
            pats = two_squares(4 * v - 2)
            if not pats: continue
            for Hs, gens in subgroups(v).items():
                if (v - 1) in Hs or len(Hs) % 2 == 0: continue
                OH = 1 + len({min(t * h % v for h in Hs) for t in range(1, v)})
                if OH > 70: continue
                for mu in mus(v, Hs):
                    tw = Twin(v, list(gens), mu)
                    for c, d in pats:
                        e = tw.estimate(c, d, rng=rng)
                        if e is None: continue
                        lE, lNA, lNC, lND = e
                        work = 10 ** lNA * 10 ** max(lNC, lND) + 10 ** min(lNC, lND)
                        if lE > -3:
                            print(f"{lE:7.2f} v={v} |H|={len(Hs)} gens={gens} mu={mu} OH={tw.stH.O} OK={tw.stK.O} classesK={tw.stK.C} "
                                  f"rows=(1,1,{c},{d}) log10N=({lNA:.1f},{lNC:.1f},{lND:.1f}) work={work:.1e}", flush=True)


def run(v, hgens, mu, c, d, cap=1e11, K=2):
    tw = Twin(v, hgens, mu); stH, stK = tw.stH, tw.stK; pm = v % 4
    NC, ND = int(stK.pool(c)), int(stK.pool(d))
    if NC > ND: c, d, NC, ND = d, c, ND, NC
    nA = 2 ** len(tw.pairs)
    work = nA * ND + NC
    print(f"v={v} H={hgens} mu={mu} OH={stH.O} OK={stK.O} classes={stK.C} rows=(1,1,{c},{d}) pools A={nA} C={NC} D={ND} work={work:.1e}", flush=True)
    if work > cap or len(tw.pairs) > 26: return []
    sg = 1 - 2 * ((np.arange(nA)[:, None] >> np.arange(len(tw.pairs))) & 1)
    YA = tw.A_vectors(sg.astype(np.float64))
    VA = np.empty((nA, stK.C), np.int16)
    for s0 in range(0, nA, 100000): VA[s0:s0 + 100000] = np.rint((tw.TA(YA[s0:s0 + 100000]) - 2 * pm) / 4)
    YC, YD = pool_exact(stK, c).astype(np.float64), pool_exact(stK, d).astype(np.float64)
    def vals(Y):
        V = np.empty((len(Y), stK.C), np.int16)
        for s0 in range(0, len(Y), 200000): V[s0:s0 + 200000] = np.rint((stK.class_values(Y[s0:s0 + 200000]) - pm) / 4)
        return V
    Vs = [VA, np.zeros((1, stK.C), np.int16), vals(YC), vals(YD)]
    fd, path = tempfile.mkstemp(suffix='.bin', dir=W.tmpdir()); os.close(fd)
    t = time.time()
    try:
        with open(path, 'wb') as f:
            f.write(struct.pack('iiiii', *[len(V) for V in Vs], stK.C)); f.write(np.full(stK.C, -pm, dtype=np.int16).tobytes())
            for V in Vs: f.write(np.ascontiguousarray(V).tobytes())
        out = subprocess.run([W.bin('find4b'), path, str(K)], capture_output=True, text=True, check=True).stdout
    finally:
        os.remove(path)
    total = int(out.strip().split('TOTAL')[-1])
    print(f"  {total} solutions [{time.time() - t:.0f}s]", flush=True)
    found = []
    for line in out.splitlines():
        if not line.startswith('SOL'): continue
        ia, _, ic, idd = map(int, line.split()[1:])
        a = YA[ia].astype(np.int64)[stH.orb]
        am = a[np.arange(v) * mu % v]
        seqs = [a, am, YC[ic].astype(np.int64)[stK.orb], YD[idd].astype(np.int64)[stK.orb]]
        tot = verify.paf_sum(seqs)
        if not (tot[0] == 4 * v and not tot[1:].any()): continue
        H = verify.goethals_seidel(*map(verify.circulant, seqs)).astype(np.float64)
        had = np.array_equal(H @ H.T, 4 * v * np.eye(4 * v)); skew = np.array_equal(H + H.T, 2 * np.eye(4 * v))
        print(f"  VERIFIED Hadamard={had} skew={skew}", flush=True)
        if had and skew:
            rec = dict(n=v, gens=list(hgens), twin_multiplier=mu, pattern=[int(x.sum()) for x in seqs], hadamard=True, skew=True,
                       seqs=[''.join('+' if x > 0 else '-' for x in s) for s in seqs])
            found.append(rec)
            with open(W.hits(), 'a') as f: f.write(json.dumps(rec) + '\n')
    return found


if __name__ == '__main__' and sys.argv[1] == 'run':
    v = int(sys.argv[2]); hg = [int(x) for x in sys.argv[3].split(',')]; mu = int(sys.argv[4])
    c, d = (int(x) for x in sys.argv[5].split(','))
    run(v, hg, mu, c, d, float(sys.argv[6]) if len(sys.argv) > 6 else 1e11)

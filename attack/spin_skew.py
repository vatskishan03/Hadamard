"""Spin-type skew Goethals-Seidel families (the structure of Djokovic's skew 4*397).

A is skew-type and invariant under K = <H, mu> (mu of order 3, -1 not in K); B is H-invariant and
C = B o mu, D = B o mu^2, so P_C(s) = P_B(mu s). The family is complementary iff

    P_A(s) + T_B(s) = 0  for s != 0,     T_B(s) = P_B(s) + P_B(mu s) + P_B(mu^2 s),

on the orbits of <K, -1>. The row sums are (1, b, b, b) with 1 + 3 b^2 = 4v.

  python3 spin_skew.py rank v ...           Gaussian estimate for every (H, mu)
  python3 spin_skew.py search v h mu [--restarts R] [--steps S] [--seed X]
"""
import argparse, json, math, os, sys, time
from math import gcd, log10, pi
import numpy as np
from pathlib import Path as _P; sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
import spence_lib as L
from gs_sweep_subgroups import subgroups
import verify
HERE = str(W.HERE)


def closure(gens, v):
    S = {1}
    while True:
        T = S | {s * g % v for s in S for g in gens}
        if T == S: return frozenset(S)
        S = T


class Spin:
    def __init__(self, v, hgens, mu):
        self.v, self.mu = v, mu
        G = L.AbGroup([v]); self.G = G
        P = lambda g: G.perm_from_matrix(np.array([[g]]))
        self.stH = L.Structure(G, [P(g) for g in hgens])
        self.stK = L.Structure(G, [P(g) for g in hgens] + [P(mu)])
        K = closure(list(hgens) + [mu], v)
        assert (v - 1) not in K, "-1 in K"
        self.H, self.K = closure(list(hgens), v), K
        stK, stH = self.stK, self.stH
        rep = np.array([np.flatnonzero(stK.orb == o)[0] for o in range(stK.O)])
        partner = stK.orb[G.neg()[rep]]
        self.pairs = [(o, partner[o]) for o in range(1, stK.O) if o < partner[o]]
        # T_B on K-class c = sum_k P_B at the H-class of mu^k s_c
        self.M = np.zeros((stK.C, stH.O, stH.O))
        for c, s in enumerate(stK.reps):
            for k in range(3):
                self.M[c] += stH.S[stH.cls[s * pow(mu, k, v) % v]]
        self.M = (self.M + self.M.transpose(0, 2, 1)) / 2          # same quadratic forms, symmetric
        self.SA = (stK.S + stK.S.transpose(0, 2, 1)) / 2
        self.b = math.isqrt((4 * v - 1) // 3)
        assert 1 + 3 * self.b ** 2 == 4 * v

    def A_vectors(self, s):
        Y = np.zeros((len(s), self.stK.O)); Y[:, 0] = 1
        for j, (o, q) in enumerate(self.pairs): Y[:, o] = s[:, j]; Y[:, q] = -s[:, j]
        return Y

    def PA(self, YA):
        return self.stK.class_values(YA)

    def TB(self, YB):
        return np.einsum("bo,cop,bp->bc", YB, self.M, YB, optimize=True)

    def estimate(self, samples=3000, rng=None):
        rng = rng or np.random.default_rng(1)
        nA = len(self.pairs)
        s = np.where(rng.random((samples, nA)) < 0.5, 1.0, -1.0)
        VA = self.PA(self.A_vectors(s))
        NB = self.stH.pool(self.b)
        if NB == 0: return None
        VB = self.TB(self.stH.sample(self.b, samples, rng))
        V = VA + VB
        lp = []
        for c in range(self.stK.C):
            var, mu = V[:, c].var(), V[:, c].mean()
            if var < 1e-9:
                if abs(mu) > 1e-6: return None
                continue
            lp.append((self.stK.csize[c] * math.sqrt(var), -0.5 * log10(2 * pi * var) - mu ** 2 / (2 * var) / math.log(10) + log10(4)))
        lp.sort()
        return nA * log10(2) + log10(NB) + sum(t for _, t in lp[:-1]), nA * log10(2), log10(NB)

    def search(self, restarts, steps, rng, verbose=True):
        """Tabu search on (A, B): A flips a negation pair, B swaps two equal-size orbits of opposite sign.
        Objective: sum over K-classes of class size * (P_A + T_B)^2."""
        stH, stK = self.stH, self.stK
        w = stK.csize.astype(float)
        SA = self.SA
        sizes = stH.sizes
        for r in range(restarts):
            B = stH.sample(self.b, 1, rng)[0]
            sA = np.where(rng.random(len(self.pairs)) < 0.5, 1.0, -1.0)
            A = self.A_vectors(sA[None])[0]
            gB = np.einsum("cop,p->co", self.M, B); TB = gB @ B
            gA = np.einsum("cop,p->co", SA, A); PA = gA @ A
            tabu = {}
            best = None
            for it in range(steps):
                R = PA + TB
                f = float(w @ R ** 2)
                if best is None or f < best: best = f
                if f == 0:
                    return A, B, r, it
                # B swaps: i (+) <-> j (-) of equal size, delta on each class
                dB = -4 * B[None, :] * gB + 4 * np.einsum("coo->co", self.M)        # (C, O) single flips
                cand = []
                plus, minus = np.flatnonzero(B > 0), np.flatnonzero(B < 0)
                plus = plus[plus > 0]; minus = minus[minus > 0]
                D = dB[:, plus][:, :, None] + dB[:, minus][:, None, :] + 8 * B[plus][None, :, None] * B[minus][None, None, :] * self.M[:, plus][:, :, minus]
                same = sizes[plus][:, None] == sizes[minus][None, :]
                fB = np.einsum("c,cij->ij", w, (R[:, None, None] + D) ** 2)
                fB[~same] = np.inf
                tabu = {k: t for k, t in tabu.items() if t > it}
                for (i, j), t in tabu.items():
                    if True:
                        ii = np.flatnonzero(plus == i); jj = np.flatnonzero(minus == j)
                        if len(ii) and len(jj): fB[ii[0], jj[0]] = np.inf
                # A pair flips
                dA = np.zeros((stK.C, len(self.pairs)))
                for k, (o, q) in enumerate(self.pairs):
                    x = A.copy(); x[o] = -x[o]; x[q] = -x[q]
                    dA[:, k] = np.einsum("cop,o,p->c", SA, x, x) - PA
                fA = np.einsum("c,ck->k", w, (R[:, None] + dA) ** 2)
                kb = np.unravel_index(np.argmin(fB), fB.shape)
                ka = int(np.argmin(fA))
                noise = rng.random() < 0.05
                if noise:
                    ii = rng.integers(len(plus)); jj = np.flatnonzero(same[ii])
                    if len(jj) == 0: continue
                    kb = (ii, rng.choice(jj)); fBmin = fB[kb]; take_B = True
                else:
                    take_B = fB[kb] <= fA[ka]
                if take_B:
                    i, j = plus[kb[0]], minus[kb[1]]
                    for x in (i, j):
                        gB -= 2 * B[x] * self.M[:, :, x]; B[x] = -B[x]
                    TB = gB @ B
                    tabu[(j, i)] = it + 7
                else:
                    o, q = self.pairs[ka]
                    for x in (o, q):
                        gA -= 2 * A[x] * SA[:, :, x]; A[x] = -A[x]
                    PA = gA @ A
            if verbose and r % 10 == 0: print(f"  restart {r}: best objective {best:.0f}", flush=True)
        return None

    def sequences(self, A, B):
        a = A.astype(np.int64)[self.stK.orb]
        b = B.astype(np.int64)[self.stH.orb]
        t = np.arange(self.v)
        c = b[t * self.mu % self.v]; d = b[t * self.mu * self.mu % self.v]
        return [a, b, c, d]


def gs_local(v, hgens, pattern, restarts, steps, seed, verbose=True):
    """Plain skew GS (no spin) by the same tabu search: B, C, D stacked into one vector with a
    block-diagonal quadratic form; swaps stay inside a block (sizes tagged by block)."""
    G = L.AbGroup([v]); P = lambda g: G.perm_from_matrix(np.array([[g]]))
    st = L.Structure(G, [P(g) for g in hgens])
    rep = np.array([np.flatnonzero(st.orb == o)[0] for o in range(st.O)])
    partner = st.orb[G.neg()[rep]]
    assert not np.any(partner[1:] == np.arange(1, st.O)), "-1 fixes an orbit"
    pairs = [(o, partner[o]) for o in range(1, st.O) if o < partner[o]]
    O, C = st.O, st.C
    S2 = np.rint(st.S + st.S.transpose(0, 2, 1)).astype(np.int32)
    M2 = np.zeros((C, 3 * O, 3 * O), np.int32)
    for k in range(3): M2[:, k * O:(k + 1) * O, k * O:(k + 1) * O] = S2
    rng = np.random.default_rng(seed)
    X0 = np.concatenate([st.sample(r, 1, rng)[0] for r in pattern]).astype(np.int32)
    tags = np.concatenate([st.sizes + 100000 * k for k in range(3)]).astype(np.int32)
    path = os.path.join(W.tmpdir(), f'gsls_{v}_{seed}.bin')
    with open(path, 'wb') as f:
        np.array([3 * O, C, O, len(pairs)], np.int32).tofile(f)
        st.csize.astype(np.int32).tofile(f); tags.tofile(f)
        M2.tofile(f); S2.tofile(f); np.array(pairs, np.int32).ravel().tofile(f); X0.tofile(f)
    import subprocess
    t = time.time()
    out = subprocess.run([W.bin(os.environ.get('SPIN_LS', 'spin_ls')), path, str(restarts), str(steps), str(seed)] + os.environ.get('SPIN_LS_ARGS', '').split(), capture_output=True, text=True)
    os.remove(path)
    if verbose: print(f"v={v} gens={hgens} rows=(1,{pattern}) orbits={O} classes={C}:", out.stderr.strip().splitlines()[-1], f"[{time.time()-t:.0f}s]", flush=True)
    sol = [l for l in out.stdout.splitlines() if l.startswith('SOL')]
    if not sol: return None
    tok = sol[0].split(); ia = tok.index('A'); ib = tok.index('B')
    A = np.array(tok[ia + 1:ib], dtype=np.int64); X = np.array(tok[ib + 1:], dtype=np.int64)
    seqs = [A[st.orb]] + [X[k * O:(k + 1) * O][st.orb] for k in range(3)]
    tot = verify.paf_sum(seqs)
    ok = tot[0] == 4 * v and not tot[1:].any()
    H = verify.goethals_seidel(*map(verify.circulant, seqs)).astype(np.float64)
    had = np.array_equal(H @ H.T, 4 * v * np.eye(4 * v)); skew = np.array_equal(H + H.T, 2 * np.eye(4 * v))
    print(f"FOUND: complementary={ok} Hadamard={had} skew={skew}", flush=True)
    if ok and had and skew:
        with open(W.hits(), 'a') as f:
            f.write(json.dumps(dict(n=v, gens=list(hgens), pattern=[int(x.sum()) for x in seqs], hadamard=True, skew=True,
                                    seqs=[''.join('+' if x > 0 else '-' for x in s) for s in seqs])) + '\n')
        return seqs
    return None


def spin_exhaustive(v, hgens, mu, maxpool=2e7):
    """Exact 2-list join for small spin structures: every skew A (K-invariant) against every B."""
    sp = Spin(v, hgens, mu)
    NB = sp.stH.pool(sp.b)
    if NB > maxpool or len(sp.pairs) > 24 or sp.stH.O > 44: return None
    import skew_ktree as K
    YB = K.pool_exact(sp.stH, sp.b).astype(np.float64)
    s = 1 - 2 * ((np.arange(2 ** len(sp.pairs))[:, None] >> np.arange(len(sp.pairs))) & 1)
    YA = sp.A_vectors(s.astype(np.float64))
    PA = np.rint(sp.PA(YA)).astype(np.int64)
    table = {}
    for i, row in enumerate(PA): table.setdefault(tuple(-row), []).append(i)
    hits = []
    for s0 in range(0, len(YB), 200000):
        TB = np.rint(sp.TB(YB[s0:s0 + 200000])).astype(np.int64)
        for j, row in enumerate(TB):
            for i in table.get(tuple(row), ()):
                hits.append((YA[i], YB[s0 + j]))
    return sp, hits


def order3(v, K0):
    """Multipliers mu of order 3 modulo K0: mu not in K0 and mu^3 in K0 (mu^3 = 1 is not required: the blocks
    B, B o mu, B o mu^2 close up because B is K0-invariant). One representative per coset mu K0, mu and mu^-1
    counted separately (they give mirror-image structures)."""
    out, seen = [], set()
    for u in range(2, v):
        if gcd(u, v) == 1 and u not in K0 and pow(u, 3, v) in K0:
            c = frozenset(u * k % v for k in K0)
            if c not in seen: seen.add(c); out.append(u)
    return out

if __name__ == '__main__':
    if sys.argv[1] == 'rank':
        rng = np.random.default_rng(1)
        for v in map(int, sys.argv[2:]):
            for Hs, gens in subgroups(v).items():
                if (v - 1) in Hs or len(Hs) % 2 == 0: continue
                for mu in order3(v, Hs):
                    K = closure(list(gens) + [mu], v)
                    if (v - 1) in K: continue
                    try:
                        sp = Spin(v, list(gens), mu)
                    except AssertionError:
                        continue
                    if sp.stH.O > 120: continue
                    e = sp.estimate(rng=rng)
                    if e is None: continue
                    print(f"v={v} |H|={len(Hs)} gens={gens} mu={mu} OH={sp.stH.O} OK={sp.stK.O} classesK={sp.stK.C} "
                          f"log10E={e[0]:.2f} log10NA={e[1]:.1f} log10NB={e[2]:.1f}", flush=True)
                    pass
    elif sys.argv[1] == 'gsearch':
        v = int(sys.argv[2]); hg = [int(x) for x in sys.argv[3].split(',')]; pat = [int(x) for x in sys.argv[4].split(',')]
        gs_local(v, hg, pat, int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]) if len(sys.argv) > 7 else 1)
    elif sys.argv[1] == 'exhaust':
        rng = np.random.default_rng(1)
        for v in map(int, sys.argv[2:]):
            for Hs, gens in subgroups(v).items():
                if (v - 1) in Hs or len(Hs) % 2 == 0: continue
                for mu in order3(v, Hs):
                    if (v - 1) in closure(list(gens) + [mu], v): continue
                    try:
                        r = spin_exhaustive(v, list(gens), mu)
                    except AssertionError:
                        continue
                    if r is None: continue
                    sp, hits = r
                    print(f"v={v} |H|={len(Hs)} gens={gens} mu={mu} A={2 ** len(sp.pairs)} B={sp.stH.pool(sp.b):.0f}: {len(hits)} solutions", flush=True)
                    if hits:
                        A, B = hits[0]; seqs = sp.sequences(A, B)
                        tot = verify.paf_sum(seqs); ok = tot[0] == 4 * v and not tot[1:].any()
                        H = verify.goethals_seidel(*map(verify.circulant, seqs)).astype(np.float64)
                        had = np.array_equal(H @ H.T, 4 * v * np.eye(4 * v)); skew = np.array_equal(H + H.T, 2 * np.eye(4 * v))
                        print(f"  VERIFIED complementary={ok} Hadamard={had} skew={skew}", flush=True)
                        if ok and had and skew:
                            with open(W.hits(), 'a') as f:
                                f.write(json.dumps(dict(n=v, gens=list(gens), spin_multiplier=mu, pattern=[int(x.sum()) for x in seqs],
                                                        hadamard=True, skew=True, seqs=[''.join('+' if x > 0 else '-' for x in q) for q in seqs])) + '\n')
    elif sys.argv[1] == 'csearch':
        ap = argparse.ArgumentParser(); ap.add_argument('cmd'); ap.add_argument('v', type=int); ap.add_argument('h')
        ap.add_argument('mu', type=int); ap.add_argument('--restarts', type=int, default=100000)
        ap.add_argument('--steps', type=int, default=2000); ap.add_argument('--seed', type=int, default=1)
        a = ap.parse_args()
        sp = Spin(a.v, [int(x) for x in a.h.split(',')], a.mu)
        stH, stK = sp.stH, sp.stK
        print(f"v={a.v} |H|={len(sp.H)} |K|={len(sp.K)} OH={stH.O} OK={stK.O} classesK={stK.C} b={sp.b}", flush=True)
        M2 = np.rint(2 * sp.M).astype(np.int32); S2 = np.rint(2 * sp.SA).astype(np.int32)
        assert np.allclose(M2, 2 * sp.M) and np.allclose(S2, 2 * sp.SA)
        B0 = stH.sample(sp.b, 1, np.random.default_rng(a.seed))[0].astype(np.int32)
        path = os.path.join(W.tmpdir(), f'spin_{a.v}_{a.seed}.bin')
        with open(path, 'wb') as f:
            np.array([stH.O, stK.C, stK.O, len(sp.pairs)], np.int32).tofile(f)
            stK.csize.astype(np.int32).tofile(f); stH.sizes.astype(np.int32).tofile(f)
            M2.tofile(f); S2.tofile(f); np.array(sp.pairs, np.int32).ravel().tofile(f); B0.tofile(f)
        import subprocess
        t = time.time()
        if os.environ.get('SPIN_WALK'):      # neighbour-lookup walk: SPIN_WALK = seconds
            out = subprocess.run([W.bin('spin_walk'), path, os.environ['SPIN_WALK'], str(a.seed)], capture_output=True, text=True)
        else:
            out = subprocess.run([W.bin(os.environ.get('SPIN_LS', 'spin_ls')), path, str(a.restarts), str(a.steps), str(a.seed)]
                                 + os.environ.get('SPIN_LS_ARGS', '').split(), capture_output=True, text=True)
        os.remove(path)
        print(out.stderr.strip().splitlines()[-1], f"[{time.time()-t:.0f}s]", flush=True)
        sol = [l for l in out.stdout.splitlines() if l.startswith('SOL')]
        if not sol: sys.exit(0)
        tok = sol[0].split(); ia = tok.index('A'); ib = tok.index('B')
        A = np.array(tok[ia + 1:ib], dtype=float); B = np.array(tok[ib + 1:], dtype=float)
        seqs = sp.sequences(A, B)
        tot = verify.paf_sum(seqs)
        ok = tot[0] == 4 * a.v and not tot[1:].any()
        H = verify.goethals_seidel(*map(verify.circulant, seqs)).astype(np.float64)
        had = np.array_equal(H @ H.T, 4 * a.v * np.eye(4 * a.v)); skew = np.array_equal(H + H.T, 2 * np.eye(4 * a.v))
        print(f"FOUND: complementary={ok} Hadamard={had} skew={skew}", flush=True)
        if ok and had and skew:
            with open(W.hits(), 'a') as f:
                f.write(json.dumps(dict(n=a.v, gens=sorted(int(x) for x in a.h.split(',')), spin_multiplier=a.mu,
                                        pattern=[int(x.sum()) for x in seqs], hadamard=True, skew=True,
                                        seqs=[''.join('+' if x > 0 else '-' for x in s) for s in seqs])) + '\n')
    else:
        ap = argparse.ArgumentParser(); ap.add_argument('cmd'); ap.add_argument('v', type=int); ap.add_argument('h')
        ap.add_argument('mu', type=int); ap.add_argument('--restarts', type=int, default=1000)
        ap.add_argument('--steps', type=int, default=3000); ap.add_argument('--seed', type=int, default=1)
        a = ap.parse_args()
        sp = Spin(a.v, [int(x) for x in a.h.split(',')], a.mu)
        print(f"v={a.v} H={sorted(sp.H)[:8]}... |H|={len(sp.H)} |K|={len(sp.K)} OH={sp.stH.O} OK={sp.stK.O} "
              f"classesK={sp.stK.C} b={sp.b}", flush=True)
        t = time.time()
        res = sp.search(a.restarts, a.steps, np.random.default_rng(a.seed))
        if res is None:
            print("no solution", flush=True); sys.exit(0)
        A, B, r, it = res
        seqs = sp.sequences(A, B)
        tot = verify.paf_sum(seqs)
        ok = tot[0] == 4 * a.v and not tot[1:].any()
        H = verify.goethals_seidel(*map(verify.circulant, seqs)).astype(np.float64)
        had = np.array_equal(H @ H.T, 4 * a.v * np.eye(4 * a.v)); skew = np.array_equal(H + H.T, 2 * np.eye(4 * a.v))
        print(f"FOUND restart {r} step {it} [{time.time()-t:.0f}s]: complementary={ok} Hadamard={had} skew={skew}", flush=True)
        if ok and had and skew:
            with open(W.hits(), 'a') as f:
                f.write(json.dumps(dict(n=a.v, gens=sorted(int(x) for x in a.h.split(',')), spin_multiplier=a.mu,
                                        pattern=[int(x.sum()) for x in seqs], hadamard=True, skew=True,
                                        seqs=[''.join('+' if x > 0 else '-' for x in s) for s in seqs])) + '\n')

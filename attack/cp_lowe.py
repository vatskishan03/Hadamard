"""Low-energy four-list k-tree for the cyclotomic products of cyclo_product.py (K_m trivial).

Every member of each non-skew block is enumerated in C (cp_enum.c). A member's key is its autocorrelation vector in the
coordinates of cyclo_product.level_vectors: level B (the compression to GF(q), classes C_0 .. C_{nb-1}) and level C
(the shifts (b, C_c), 1 <= b <= m/2). With mu_k and Sigma_k the mean and covariance of the keys of block k (from
uniform samples; block 0 exactly when it is enumerated in full), a member is kept when its key is close to the
conditional mean c_k = mu_k + Sigma_k Sigma^+ (t - sum_j mu_j) given that the four keys sum to the target t:
|W_k (key - c_k)|^2 <= tau_k with W_k = Sigma_k^(-1/2), tau_k chosen so that about L members are kept. The lists are
then joined by kt4 (attack/kt4.c), and every solution is rebuilt as four +-1 sequences and checked.

  python3 attack/cp_lowe.py v --q Q --e E [--plain] [--config C] [--L 2e7] [--keepA L] [--S 0,1,2,3,4] [--nsig N]
                            [--samples N] [--threads T] [--estimate-only] [--seed S]

The Gaussian estimate of the number of solutions inside the kept lists is printed before the join. Hits are appended
to $HADAMARD_WORK/skew_hits.jsonl (had_hits.jsonl for plain structures)."""
import argparse, itertools, json, math, os, struct, subprocess, sys, tempfile, time
import numpy as np
from math import comb
from pathlib import Path as _P
sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
import cyclo_product as CP
import verify


def rotate(e, pid, t):
    """Pattern with its classes rotated: new x_i = old x_(i - t); the value f at (a, 0) is unchanged."""
    t %= e
    if not t: return pid
    mask = pid & ((1 << e) - 1)
    return (pid & (1 << e)) | (((mask << t) | (mask >> (e - t))) & ((1 << e) - 1))


class Orbits:
    """Positions of Z_m under K_m = <kgen>, a = kgen^j(a) * rep(a). With a twist t, the multiplier acts as
    (a, y) -> (kgen a, g^t y), so the pattern at a is the pattern at rep(a) with its classes rotated by j(a) t; a
    pattern at a representative must be invariant under the rotation by (orbit length) t."""

    def __init__(self, S, twist=0):
        m, e = S.m, S.e
        self.S, self.tw = S, twist % e
        self.reps = list(S.kreps); self.ri = {r: i for i, r in enumerate(self.reps)}
        self.rep_of = [S.kreps[S.korb[a]] for a in range(m)]
        self.j = [-1] * m; self.length = {}
        for r in self.reps:
            a, j = r, 0
            while True:
                if self.j[a] < 0: self.j[a] = j
                a = a * S.kgen % m; j += 1
                if a == r: break
            self.length[r] = j
        assert all(x >= 0 for x in self.j)

    def ok(self, pid, r):
        return rotate(self.S.e, pid, self.length[r] * self.tw) == pid

    def at(self, pid, a):
        """Pattern at position a when its representative carries pid."""
        return rotate(self.S.e, pid, self.j[a] * self.tw)

    def expand(self, rep_ids):
        """Pattern ids [N x m] at all positions from the ids [N x R] at the representatives."""
        m = self.S.m
        out = np.empty((len(rep_ids), m), dtype=np.int64)
        for a in range(m):
            col = rep_ids[:, self.ri[self.rep_of[a]]]
            out[:, a] = col if self.j[a] == 0 or self.tw == 0 else np.array([self.at(int(p), a) for p in col])
        return out

    def skew_members(self, y):
        """All members of the skew block (pattern ids at every position): f = 1 and x_(j + e/2) = -x_j at 0 (and
        invariant under the twist), and the pattern at -a the partner of the pattern at a."""
        S, m = self.S, self.S.m
        orb = S.korb
        p0s = [p for p in S.skew0() if rotate(S.e, p, self.tw) == p]
        pairs, seen = [], {orb[0]}
        for r in self.reps:
            if orb[r] in seen: continue
            r2 = self.rep_of[(-r) % m]; j0 = self.j[(-r) % m]
            pairs.append((r, r2, j0)); seen |= {orb[r], orb[(-r) % m]}
        choices = [[p for p in S.patterns(*CP.fs_of(S, y[r])) if self.ok(p, r)] for r, _, _ in pairs]
        rows = []
        for p0 in p0s:
            for combo in itertools.product(*choices):
                rep = {0: p0}
                for (r, r2, j0), p in zip(pairs, combo):
                    rep[r] = p; rep[r2] = rotate(S.e, S.partner(p), -j0 * self.tw)
                rows.append([rep[r] for r in self.reps])
        return self.expand(np.array(rows, dtype=np.int64)) if rows else np.zeros((0, m), np.int64)


class Block:
    """One non-skew block: the allowed slice patterns at the orbit representatives and the pair tables of the key."""

    def __init__(self, S, y, CC, idx, orbits):
        self.S, self.y, self.O = S, y, orbits
        R = orbits.reps
        self.pats = [np.array([p for p in S.patterns(*CP.fs_of(S, y[r])) if orbits.ok(p, r)], dtype=np.int64) for r in R]
        self.n = [len(p) for p in self.pats]
        self.D = S.nb + (S.m // 2) * S.e
        self.T = {(i1, i2): np.zeros((self.n[i1], self.n[i2], self.D), dtype=np.int64) for i1 in range(len(R)) for i2 in range(i1 + 1)}
        m = S.m
        rot = {}                                              # (rep index, position) -> rows of CC for the rotated patterns
        for a in range(m):
            i = orbits.ri[orbits.rep_of[a]]
            rot[a] = idx[np.array([orbits.at(int(p), a) for p in self.pats[i]], dtype=np.int64)]
        for a1 in range(m):
            for a2 in range(a1 + 1):
                i1, i2 = orbits.ri[orbits.rep_of[a1]], orbits.ri[orbits.rep_of[a2]]
                if i1 < i2: a1_, a2_, i1_, i2_ = a2, a1, i2, i1
                else: a1_, a2_, i1_, i2_ = a1, a2, i1, i2
                T = self.T[i1_, i2_]
                X, Y = rot[a1_], rot[a2_]
                if i1_ == i2_:                                # same representative: one pattern index for both positions
                    d = np.arange(len(X)); add = np.zeros((len(X), self.D), dtype=np.int64)
                    pair_terms(S, CC, X, Y, a1_, a2_, add, diag=True)
                    T[d, d] += add
                else:
                    add = np.zeros((len(X), len(Y), self.D), dtype=np.int64)
                    pair_terms(S, CC, X, Y, a1_, a2_, add, diag=False)
                    T += add
        self.size = math.prod(self.n)

    def keys(self, I):
        """Keys of the members with pattern indices I [N x R]."""
        K = np.zeros((len(I), self.D), dtype=np.int64)
        for (r1, r2), T in self.T.items(): K += T[I[:, r1], I[:, r2]]
        return K

    def sample(self, N, rng):
        return np.stack([rng.integers(0, n, size=N) for n in self.n], axis=1)

    def ids(self, I):
        """Pattern ids [N x m] at all positions of the members with pattern indices I."""
        rep_ids = np.stack([self.pats[i][I[:, i]] for i in range(len(self.n))], axis=1)
        return self.O.expand(rep_ids)


class PrimeBlock:
    """A non-skew block over GF(q) itself (m = 1): f at 0 and one sign per cyclotomic class, with a fixed class sum s.
    The classes are split into G groups that play the part of positions: a member is one sub-pattern per group, its
    key PAF(d_c), c < nb, is a quadratic form in the class signs (cyclotomic numbers) and splits into pair terms between
    groups; the class sum is imposed through the tags of cp_enum."""

    def __init__(self, S, y, G, N):
        e, nb = S.e, S.nb
        self.S, self.f, self.s = S, *CP.fs_of(S, y[0])
        self.groups = [list(range(g * e // G, (g + 1) * e // G)) for g in range(G)]
        self.n = [1 << len(g) for g in self.groups]
        self.D = nb
        self.sub = [np.array([[1 - 2 * ((i >> b) & 1) for b in range(len(g))] for i in range(1 << len(g))], dtype=np.int64)
                    for g in self.groups]                          # sub-pattern index -> signs (bit set = -1)
        self.tags = [x.sum(1) for x in self.sub]
        h = S.half
        self.T = {}
        for a in range(G):
            for b in range(a + 1):
                ga, gb = self.groups[a], self.groups[b]
                if a == b:
                    M = N[:, ga][:, :, ga]                          # [c, i, j]
                    quad = np.einsum('pi,cij,pj->pc', self.sub[a], M, self.sub[a])
                    lin = np.zeros((self.n[a], nb), dtype=np.int64)
                    for c in range(nb):                             # f (x_c + x_(c+h)) for the classes in this group
                        for cc in (c, (c + h) % e):
                            if cc in ga: lin[:, c] += self.f * self.sub[a][:, ga.index(cc)]
                    T = np.zeros((self.n[a], self.n[a], nb), dtype=np.int64)
                    d = np.arange(self.n[a]); T[d, d] = quad + lin
                else:
                    M = N[:, ga][:, :, gb] + N[:, gb][:, :, ga].transpose(0, 2, 1)
                    T = np.einsum('pi,cij,qj->pqc', self.sub[a], M, self.sub[b])
                self.T[a, b] = T
        self.size = comb(e, (e + self.s) // 2)

    def keys(self, I):
        K = np.zeros((len(I), self.D), dtype=np.int64)
        for (r1, r2), T in self.T.items(): K += T[I[:, r1], I[:, r2]]
        return K

    def sample(self, N_, rng):
        """Uniform members: random sign vectors with class sum s, cut into groups."""
        e = self.S.e; npos = (e + self.s) // 2
        X = -np.ones((N_, e), dtype=np.int64)
        idx = np.argsort(rng.random((N_, e)), axis=1)[:, :npos]
        np.put_along_axis(X, idx, 1, axis=1)
        return np.stack([((X[:, g] < 0) * (1 << np.arange(len(g)))).sum(1) for g in self.groups], axis=1)

    def ids(self, I):
        e = self.S.e; fb = 0 if self.f == 1 else 1
        mask = np.zeros(len(I), dtype=np.int64)
        for gi, g in enumerate(self.groups):
            for b, c in enumerate(g): mask |= ((I[:, gi] >> b) & 1) << c
        return ((fb << e) | mask)[:, None]


def cyclotomic_numbers(S):
    """N[c, i, j] = #{y in C_i : y + d_c in C_j} for the class representatives d_c, c < nb."""
    q, e = S.q, S.e
    N = np.zeros((S.nb, e, e), dtype=np.int64)
    if hasattr(S, 'add'):                                          # prime-power field: indices are coordinates
        y = np.flatnonzero(S.cls >= 0)
        for c in range(S.nb):
            z = S.add[y, S.reps[1 + c]]; ok = S.cls[z] >= 0
            np.add.at(N[c], (S.cls[y[ok]], S.cls[z[ok]]), 1)
        return N
    y = np.arange(1, q)
    for c in range(S.nb):
        d = S.reps[1 + c]; z = (y + d) % q; ok = z != 0
        np.add.at(N[c], (S.cls[y[ok]], S.cls[z[ok]]), 1)
    return N


def paf_vector(S, seqs):
    """Sum of the PAFs of the sequences on the group of S (cyclic Z_v, or Z_m x Z_p^k for a prime-power fibre)."""
    return S.paf_sum(seqs) if hasattr(S, 'paf_sum') else verify.paf_sum(seqs)


def check_hit(S, seqs, gens):
    """(Goethals-Seidel condition, skew first block, invariance under the generators of H)."""
    tot = paf_vector(S, seqs)
    ok = tot[0] == 4 * S.v and not tot[1:].any()
    skew_ok = True
    if S.skew:
        neg = S.neg() if hasattr(S, 'neg') else (-np.arange(S.v)) % S.v
        x = seqs[0]; skew_ok = x[0] == 1 and bool(np.all(x[neg[1:]] == -x[1:]))
    inv = True
    for g_ in gens:
        perm = S.perm_of(g_) if hasattr(S, 'perm_of') else (np.arange(S.v) * g_) % S.v
        inv = inv and all(np.array_equal(x[perm], x) for x in seqs)
    return bool(ok), bool(skew_ok), bool(inv)


def prime_key(S, N, X, f):
    """Keys PAF(d_c), c < nb, of sign vectors X [n x e] with value f at 0."""
    h = S.half
    K = np.einsum('ni,cij,nj->nc', X, N, X)
    for c in range(S.nb): K[:, c] += f * (X[:, c] + X[:, (c + h) % S.e])
    return K


def prepare_prime(S, y, samples, rng, G=4):
    """prepare() for m = 1: the skew block (2^(e/2) members) is enumerated directly, the other blocks are PrimeBlocks."""
    e, h = S.e, S.e // 2
    N = cyclotomic_numbers(S)
    blocks, mus, sigs, full, Pm0 = [], [], [], [], None
    for k in range(4):
        if S.skew and k == 0:
            half = np.array([[1 - 2 * ((i >> b) & 1) for b in range(h)] for i in range(1 << h)], dtype=np.int64)
            X = np.concatenate([half, -half], axis=1)
            K0 = np.concatenate([prime_key(S, N, X[i:i + 65536], 1) for i in range(0, len(X), 65536)])
            Pm0 = (((X < 0) * (1 << np.arange(e))).sum(1))[:, None]
            blocks.append(None); full.append(K0); mus.append(K0.mean(0)); sigs.append(np.cov(K0.T.astype(float)))
            continue
        blk = PrimeBlock(S, y[k], G, N)
        I = blk.sample(samples, rng); K = blk.keys(I)
        if k == 1 or (not S.skew and k == 0):                       # cross-check against the PAF of the sequences
            for t in range(5):
                x = S.sequence(blk.ids(I[t:t + 1])[0])
                paf = paf_vector(S, [x])
                assert np.array_equal(paf[[S.reps[1 + c] for c in range(S.nb)]], K[t]), 'prime pair tables disagree with the PAF'
        blocks.append((blk, I, K)); full.append(None); mus.append(K.mean(0)); sigs.append(np.cov(K.T.astype(float)))
    D = S.nb
    sizes = [len(f) if f is not None else b[0].size for f, b in zip(full, blocks)]
    est0, rank0, _ = gauss_estimate([f if f is not None else b[2] for f, b in zip(full, blocks)], np.zeros(D))
    est0 += sum(math.log10(b[0].size / len(b[2])) for b in blocks if b is not None)
    return dict(blocks=blocks, full=full, mus=mus, sigs=sigs, Pm0=Pm0, D=D, sizes=sizes, est0=est0, rank0=rank0)


def stage1_prime(S):
    """Stage-1 classes for m = 1: compressions y = f + F s of the blocks (skew block y = 1) with sum of squares 4v, up to
    negation and order of the non-skew blocks."""
    e, F, v = S.e, S.F, S.v
    vals = sorted({abs(f + F * s) for s in range(-e, e + 1, 2) for f in (1, -1)})
    lead = [1] if S.skew else vals
    out = []
    for y0 in lead:
        for i, a in enumerate(vals):
            for j, b in enumerate(vals[i:], i):
                c2 = 4 * v - y0 * y0 - a * a - b * b
                if c2 < b * b: continue
                c = math.isqrt(c2)
                if c * c != c2 or c not in vals: continue
                if not S.skew and not (y0 <= a): continue
                out.append(((y0,), (a,), (b,), (c,)))
    return out


def pair_terms(S, CC, X, Y, a1, a2, out, diag):
    """Add to out the key contribution of the positions a1, a2 carrying the CC rows X (at a1) and Y (at a2): every
    pair (X[i], Y[k]) into out[i, k], or only (X[i], Y[i]) into out[i] when diag."""
    m, e, nb = S.m, S.e, S.nb
    if diag: A = CC[X, Y].astype(np.int64); B = CC[Y, X].astype(np.int64)
    else: A = CC[X][:, Y].astype(np.int64); B = CC[Y][:, X].astype(np.int64).swapaxes(0, 1)
    if a1 == a2:
        out[..., :nb] += A[..., 1:1 + nb]
        return
    out[..., :nb] += A[..., 1:1 + nb] + B[..., 1:1 + nb]
    b1, b2 = (a1 - a2) % m, (a2 - a1) % m
    if 1 <= b1 <= m // 2: out[..., nb + (b1 - 1) * e: nb + b1 * e] += A[..., 1:1 + e]
    if 1 <= b2 <= m // 2: out[..., nb + (b2 - 1) * e: nb + b2 * e] += B[..., 1:1 + e]


def whitening(Sig, tol=1e-9):
    ev, U = np.linalg.eigh(Sig)
    keep = ev > tol * ev.max()
    return (U[:, keep] / np.sqrt(ev[keep])).T, int(keep.sum())


def gauss_estimate(lists, t, cap=10 ** 6, rng=None, sizes=None):
    """log10 of the expected number of quadruples summing to t, from the means and covariances of the four key lists
    (lattice 4 per coordinate of the range). With sizes, the lists are samples of lists of those sizes."""
    rng = rng or np.random.default_rng(1)
    mus, sigs, logN = [], [], 0.0
    for i, K in enumerate(lists):
        X = K if len(K) <= cap else K[rng.choice(len(K), cap, replace=False)]
        X = X.astype(float); mus.append(X.mean(0)); sigs.append(np.cov(X.T))
        logN += math.log10(sizes[i] if sizes is not None else len(K))
    mu = sum(mus); Sig = sum(sigs)
    ev, U = np.linalg.eigh(Sig); keep = ev > 1e-9 * ev.max()
    z = U[:, keep].T @ (np.asarray(t, float) - mu)
    off = float(np.abs(U[:, ~keep].T @ (np.asarray(t, float) - mu)).max()) if (~keep).any() else 0.0
    dens = sum(-0.5 * math.log10(2 * math.pi * x) for x in ev[keep]) - 0.5 * float((z * z / ev[keep]).sum()) / math.log(10)
    return logN + dens + int(keep.sum()) * math.log10(4), int(keep.sum()), off


def run_enum(blk, Wk, ck, tau, path_out, threads=None, stream=False):
    """Write the cp_enum input for block blk and run it; returns the kept records as (masks, keys, bits, message)
    (with stream, the records stay in path_out and masks, keys are None)."""
    m, D = len(blk.n), blk.D
    Z = Wk.shape[0]
    bits = max(1, max(int(n - 1).bit_length() for n in blk.n))
    fd, path_in = tempfile.mkstemp(suffix='.cpe', dir=W.tmpdir()); os.close(fd)
    try:
        with open(path_in, 'wb') as f:
            f.write(struct.pack('4i', m, D, Z, bits)); f.write(np.array(blk.n, np.int32).tobytes())
            f.write(struct.pack('f', tau)); f.write((-(Wk @ ck)).astype(np.float32).tobytes())
            f.write(np.zeros(D, np.int32).tobytes())
            for r1 in range(m):
                for r2 in range(r1 + 1):
                    T = blk.T[r1, r2]
                    f.write(T.astype(np.int32).tobytes())
                    f.write((T.astype(np.float64) @ Wk.T).astype(np.float32).tobytes())
            if isinstance(blk, PrimeBlock):                        # class sum s through the tags
                f.write(struct.pack('i', blk.s))
                for tg in blk.tags: f.write(np.asarray(tg, np.int32).tobytes())
        env = dict(os.environ)
        if threads: env['OMP_NUM_THREADS'] = str(threads)
        r = subprocess.run([W.bin('cp_enum'), path_in, path_out], capture_output=True, text=True, env=env, check=True)
        if stream: return None, None, bits, r.stdout.strip()
        arr = np.fromfile(path_out, dtype=np.dtype([('mask', '<u8'), ('key', '<i2', (D,))]))
        masks = arr['mask'].copy(); keys = arr['key'].copy(); del arr
        return masks, keys, bits, r.stdout.strip()
    finally:
        os.remove(path_in)


def decode(blk, mask, bits):
    """Pattern ids at the positions 0 .. m-1 of the member with the given cp_enum mask."""
    I = np.array([[(int(mask) >> (bits * i)) & ((1 << bits) - 1) for i in range(len(blk.n))]], dtype=np.int64)
    return blk.ids(I)[0]


def colvar(L, chunk=1 << 20):
    """Column variances of an integer array, in chunks (lists of 10^7 keys do not fit in float64 at once)."""
    s = np.zeros(L.shape[1]); s2 = np.zeros(L.shape[1])
    for i in range(0, len(L), chunk):
        X = L[i:i + chunk].astype(float); s += X.sum(0); s2 += (X * X).sum(0)
    m = s / len(L)
    return s2 / len(L) - m * m


def choose_sigmas(lists, S_idx, t, nsig, M, rng, sizes=None):
    """Split values sigma of the merge classes, ranked by the expected number of solutions per unit of kt4 work.
    AB(sigma) and CD(sigma) (pairs (a, b) with (a + b)_S = sigma and pairs (c, d) with (c + d)_S = t_S - sigma) are
    estimated from M random pairs per side, which respects sublattices of the key values. With X = a + b and Y = c + d,
    the solutions with split sigma number AB CD rho(sigma), rho = N(sigma; m, C) / (N_X(sigma) N_Y(t_S - sigma)), where
    N(m, C) is the Gaussian law of X_S given X + Y = t: the other classes are correlated with the merge classes, so
    the solutions do not sit where AB CD is largest. kt4 spends about AB + CD + |A| + |C| operations per value."""
    S_idx = list(S_idx)
    def stats(L):
        X = L if len(L) <= 400000 else L[rng.choice(len(L), 400000, replace=False)]
        X = X.astype(float); return X.mean(0), np.cov(X.T)
    st = [stats(L) for L in lists]
    muX, SX = st[0][0] + st[1][0], st[0][1] + st[1][1]; muY, SY = st[2][0] + st[3][0], st[2][1] + st[3][1]
    G = SX @ np.linalg.pinv(SX + SY, rcond=1e-10)
    m = (muX + G @ (t - muX - muY))[S_idx]; C = (SX - G @ SX)[np.ix_(S_idx, S_idx)]
    def logN(Z, mean, cov):
        ev, U = np.linalg.eigh(cov); keep = ev > 1e-9 * ev.max()
        w = (Z - mean) @ U[:, keep]
        return -0.5 * (w * w / ev[keep]).sum(1) - 0.5 * np.log(ev[keep]).sum()
    def sums(L1, L2, sign):
        i, j = rng.integers(0, len(L1), M), rng.integers(0, len(L2), M)
        X = L1[i][:, S_idx].astype(np.int64) + L2[j][:, S_idx].astype(np.int64)
        return X if sign > 0 else t[S_idx][None, :] - X
    ab = sums(lists[0], lists[1], 1); cd = sums(lists[2], lists[3], -1)
    lo = np.minimum(ab.min(0), cd.min(0)); span = np.maximum(ab.max(0), cd.max(0)) - lo + 1
    enc = lambda X: np.ravel_multi_index((X - lo).T, span)
    ua, ca = np.unique(enc(ab), return_counts=True); uc, cc = np.unique(enc(cd), return_counts=True)
    common, ia, ic = np.intersect1d(ua, uc, return_indices=True)
    n_ = sizes if sizes is not None else [len(l) for l in lists]
    nab = ca[ia] / M * n_[0] * n_[1]; ncd = cc[ic] / M * n_[2] * n_[3]
    sig = (np.array(np.unravel_index(common, span)).T + lo).astype(float)
    lrho = logN(sig, m, C) - logN(sig, muX[S_idx], SX[np.ix_(S_idx, S_idx)]) - logN(t[S_idx] - sig, muY[S_idx], SY[np.ix_(S_idx, S_idx)])
    gain = nab * ncd * np.exp(lrho - lrho.max()); cost = nab + ncd + n_[0] + n_[2]
    order = np.argsort(-gain / cost)[:nsig]
    share = gain[order].sum() / gain.sum() if gain.sum() else 0.0
    return sig[order].astype(np.int64), dict(candidates=len(common), share=float(share), ab=float(nab[order].sum()), cd=float(ncd[order].sum()))


def skip_done(a, sigmas, info, log):
    """Drop the split values listed in a.skip_sigmas (one tuple per line, as kt4 logs them: 'kt4: sigma i/N (s0,..)' or
    plain integers) and keep at most a.nsig of the rest; used to resume a run whose lists are rebuilt identically."""
    if not a.skip_sigmas: return sigmas[:a.nsig]
    import re
    done = set()
    for line in open(a.skip_sigmas):
        m = re.search(r'\(([-0-9,]+)\)', line)
        vals = m.group(1).split(',') if m else line.split()
        if vals: done.add(tuple(int(x) for x in vals))
    keep = [i for i, sg in enumerate(sigmas.tolist()) if tuple(sg) not in done]
    log(f'skipping {len(sigmas) - len(keep)} split values done before ({len(done)} listed)')
    return sigmas[keep][:a.nsig]


def stage1_cached(S, plain):
    """Stage-1 classes of the structure, cached in $HADAMARD_WORK."""
    cache = W.WORK / f'cp_stage1_{S.v}_{S.q}_{S.e}_{S.kgen}_{"plain" if plain else "skew"}.json'
    if cache.exists():
        return [tuple(tuple(y_) for y_ in qd) for qd in json.load(open(cache))]
    quads = CP.stage1(S, maxpairs=3e8, maxopts=3e6)
    json.dump(quads, open(cache, 'w'))
    return quads


def admissible(S, O, Y):
    """Rows of Y (compressions) for which every orbit representative has a slice pattern invariant under the twist."""
    if not O.tw: return np.ones(len(Y), dtype=bool)
    memo = {}
    def ok(t, r):
        key = (int(t), O.length[r])
        if key not in memo: memo[key] = any(O.ok(p, r) for p in S.patterns(*CP.fs_of(S, int(t))))
        return memo[key]
    return np.array([all(ok(y[r], r) for r in O.reps) for y in Y], dtype=bool)


def random_options(S, skewblock, M, rng, O=None):
    """M random compressions of one block (as stage1_options, without listing them all): a random set of spiked orbit
    units within the budget, random spike types and random signs elsewhere. Returns (Y, C, P) like stage1_options."""
    m, e, F = S.m, S.e, S.F
    T = 4 * m * e
    cost = lambda s_, f_: F * s_ * s_ + 2 * f_ * s_
    spikes = [(f_, s_) for s_ in range(-e, e + 1, 2) if s_ for f_ in (1, -1) if cost(s_, f_) <= T]
    orbits = [[a for a in range(m) if S.korb[a] == o] for o in range(len(S.kreps))]
    if skewblock:
        free, seen = [], {S.korb[0]}
        for o, orb in enumerate(orbits):
            if o in seen: continue
            neg = S.korb[(-orb[0]) % m]; free.append((orb, orbits[neg])); seen |= {o, neg}
        units = [2 * len(orb) for orb, _ in free]
    else:
        free = [(orb, None) for orb in orbits]; units = [len(orb) for orb, _ in free]
    R = len(free)
    SF = np.array([f_ for f_, _ in spikes]); SS = np.array([s_ for _, s_ in spikes]); SC = np.array([cost(s_, f_) for f_, s_ in spikes])
    U = np.array(units)
    nsp = rng.integers(0, min(R, 8) + 1, size=M)
    order = np.argsort(rng.random((M, R)), axis=1)
    rank = np.empty_like(order); rank[np.arange(M)[:, None], order] = np.arange(R)[None, :]
    spiked = rank < nsp[:, None]
    stype = rng.integers(0, len(spikes), size=(M, R))
    fval = np.where(spiked, SF[stype], np.where(rng.random((M, R)) < 0.5, 1, -1))
    sval = np.where(spiked, SS[stype], 0)
    C = (spiked * U[None, :] * SC[stype]).sum(1)
    yunit = fval + F * sval
    Y = np.zeros((M, m), dtype=np.int64)
    for u, (orb, neg) in enumerate(free):
        Y[:, orb] = yunit[:, u][:, None]
        if neg is not None: Y[:, neg] = -yunit[:, u][:, None]
    if skewblock: Y[:, 0] = 1
    keep = C <= T; Y, C = Y[keep], C[keep]
    ok = admissible if (O is not None and O.tw) else None
    if ok is not None:
        keep = admissible(S, O, Y); Y, C = Y[keep], C[keep]
    P = np.stack([(Y * np.roll(Y, -b, axis=1)).sum(1) for b in range(1, m // 2 + 1)], axis=1) if m > 1 else np.zeros((len(Y), 0), np.int64)
    return Y, C, P


def stage1_sample(S, n, rng, tries=2000000, O=None, maxopts=4e6, M=2000000):
    """Up to n stage-1 quadruples found by sampling: block 3's options are hashed by (cost, autocorrelations), and random
    options of blocks 0, 1, 2 from cost buckets that can still meet the budget are completed through the hash. For
    structures whose stage 1 is too large to list (every class tends to have a similar number of solutions). With a
    twist, only compressions that admit twist-invariant patterns at every representative are used."""
    T = 4 * S.m * S.e
    try:
        L = [CP.stage1_options(S, S.skew and k == 0, maxopts) for k in range(4)]
        if O is not None and O.tw:
            L = [tuple(x[admissible(S, O, Lk[0])] for x in Lk) for Lk in L]
    except RuntimeError:                                          # too many to list: random options instead
        L = [random_options(S, S.skew and k == 0, M if k == 3 else M // 4, rng, O) for k in range(4)]
    Y3, C3, P3 = L[3]
    table = {}
    for i in range(len(Y3)): table.setdefault((int(C3[i]), P3[i].tobytes()), []).append(i)
    byc = [{int(c): np.flatnonzero(Lk[1] == c) for c in np.unique(Lk[1])} for Lk in L]
    combos = [(c0, c1, c2, T - c0 - c1 - c2) for c0 in byc[0] for c1 in byc[1] for c2 in byc[2] if T - c0 - c1 - c2 in byc[3]]
    if not combos: return []
    w = np.array([len(byc[0][a]) * len(byc[1][b]) * len(byc[2][c]) * len(byc[3][d]) for a, b, c, d in combos], dtype=float)
    w /= w.sum()
    found, seen = [], set()
    for _ in range(tries):
        c0, c1, c2, c3 = combos[rng.choice(len(combos), p=w)]
        a = rng.choice(byc[0][c0]); b = rng.choice(byc[1][c1]); c = rng.choice(byc[2][c2])
        need = -(L[0][2][a] + L[1][2][b] + L[2][2][c])
        hits = table.get((c3, need.astype(np.int64).tobytes()))
        if not hits: continue
        d = hits[rng.integers(len(hits))]
        qd = tuple(tuple(int(x) for x in L[k][0][j]) for k, j in enumerate((a, b, c, d)))
        if qd in seen: continue
        seen.add(qd); found.append(qd)
        if len(found) >= n: break
    return found


def generators(S, O):
    """Generators of the multiplier group H: h (e-th powers on GF(q)) and the (twisted) generator of K_m. For a
    prime-power fibre (cp_field.FieldStructure) they are coordinate matrices of automorphisms of Z_m x Z_p^k."""
    if hasattr(S, 'automorphisms'): return S.automorphisms(O.tw)
    if len(S.K) == 1: return [S.h]
    return [S.h, next(t for t in range(1, S.v) if t % S.m == S.kgen and t % S.q == pow(S.g, O.tw, S.q))]


def prepare(S, O, y, samples, rng):
    """Pattern table, blocks and key statistics of the stage-1 class y (block 0 enumerated in full when skew)."""
    pats = set()
    for k in range(4):
        if S.skew and k == 0: continue
        for t in range(S.m):
            r = O.rep_of[t]
            pats |= {O.at(p, t) for p in S.patterns(*CP.fs_of(S, y[k][r])) if O.ok(p, r)}
    Pm0 = O.skew_members(y[0]) if S.skew else None
    if Pm0 is not None: pats |= set(np.unique(Pm0).tolist())
    used = sorted(pats)
    idx = np.full(1 << (S.e + 1), -1, dtype=np.int64); idx[used] = np.arange(len(used))
    CC = S.cc_table(used)
    S.bkey = []; S.dkey = S.nb
    blocks, mus, sigs, full = [], [], [], []
    for k in range(4):
        if S.skew and k == 0:
            B_, C_ = CP.level_vectors(S, Pm0.astype(np.int64), idx, CC)
            K0 = np.concatenate([B_, C_], axis=1).astype(np.int64)
            blocks.append(None); full.append(K0)
            mus.append(K0.mean(0)); sigs.append(np.cov(K0.T.astype(float)))
            continue
        blk = Block(S, y[k], CC, idx, O)
        I = blk.sample(samples, rng); K = blk.keys(I)
        if k == 1 or (not S.skew and k == 0):                       # cross-check the pair tables against level_vectors
            Bv, Cv = CP.level_vectors(S, blk.ids(I[:2000]), idx, CC)
            assert np.array_equal(np.concatenate([Bv, Cv], axis=1), K[:2000]), 'pair tables disagree with level_vectors'
        blocks.append((blk, I, K)); full.append(None)
        mus.append(K.mean(0)); sigs.append(np.cov(K.T.astype(float)))
    D = len(mus[0])
    sizes = [len(f) if f is not None else b[0].size for f, b in zip(full, blocks)]
    est0, rank0, _ = gauss_estimate([f if f is not None else b[2] for f, b in zip(full, blocks)], np.zeros(D))
    # rescale the sampled lists to the full ones
    est0 += sum(math.log10(b[0].size / len(b[2])) for b in blocks if b is not None)
    return dict(blocks=blocks, full=full, mus=mus, sigs=sigs, Pm0=Pm0, D=D, sizes=sizes, est0=est0, rank0=rank0, idx=idx, CC=CC)


def estimate_in_lists(P, Lk):
    """log10 of the Gaussian estimate of the solutions inside low-energy lists of about Lk(k) members per non-skew
    block (from the sampled members below the energy thresholds)."""
    blocks, full, mus, sigs, D = P['blocks'], P['full'], P['mus'], P['sigs'], P['D']
    t = np.zeros(D); Sp = np.linalg.pinv(sum(sigs), rcond=1e-10); mu = sum(mus)
    lists, corr = [], 0.0
    for k in range(4):
        ck = mus[k] + sigs[k] @ Sp @ (t - mu); Wk, _ = whitening(sigs[k])
        if blocks[k] is None: lists.append(full[k]); continue
        blk, I, K = blocks[k]
        en = (((K - ck) @ Wk.T) ** 2).sum(1)
        frac = min(1.0, Lk(k) / blk.size)
        sel = en <= (float(np.quantile(en, frac)) if frac < 1 else float(en.max() + 1))
        lists.append(K[sel]); corr += math.log10(min(Lk(k), blk.size) / sel.sum())
    est, rank, off = gauss_estimate(lists, t)
    return est + corr, rank, off


def stream_search(a, S, P, gens, t, Sp, mu, Lk, est0, log, rng, nsamp=2000000, chunk=1 << 22):
    """The search of main() with the lists streamed to the kt4 file: in memory stay the masks, column statistics
    and a random sample of each list (for the estimate and the split values); kt4 maps the file."""
    blocks, full, mus, sigs, Pm0, D = P['blocks'], P['full'], P['mus'], P['sigs'], P['Pm0'], P['D']
    srng = np.random.default_rng(12345)
    fd, path = tempfile.mkstemp(suffix='.kt4', dir=W.tmpdir()); os.close(fd)
    sizes, samples, recs, varsum, absmax = [], [], [], np.zeros(D), np.zeros(D, np.int64)
    found = []
    try:
        with open(path, 'wb') as fk:
            fk.write(struct.pack('5i', 0, 0, 0, 0, D)); fk.write(np.rint(t).astype(np.int16).tobytes())
            for k in range(4):
                ck = mus[k] + sigs[k] @ Sp @ (t - mu)
                Wk, rk = whitening(sigs[k])
                if blocks[k] is None:
                    K0 = full[k]
                    keep = np.arange(len(K0))
                    if a.keepA and a.keepA < len(K0):
                        en = (((K0 - ck) @ Wk.T) ** 2).sum(1); keep = np.argsort(en)[:int(a.keepA)]
                    L0 = K0[keep]
                    fk.write(np.ascontiguousarray(L0, dtype=np.int16).tobytes())
                    sizes.append(len(L0)); recs.append(('full', keep))
                    samples.append(L0 if len(L0) <= nsamp else L0[srng.choice(len(L0), nsamp, replace=False)])
                    varsum += colvar(L0); absmax = np.maximum(absmax, np.abs(L0).max(0))
                    log(f'block {k}: {len(L0)} of {len(K0)} members (skew block)')
                    continue
                blk, I, K = blocks[k]
                en = (((K - ck) @ Wk.T) ** 2).sum(1)
                frac = min(1.0, Lk(k) / blk.size)
                tau = float(np.quantile(en, frac)) if frac < 1 else 3.0e38
                fd2, pout = tempfile.mkstemp(suffix=f'.b{k}', dir=W.tmpdir()); os.close(fd2)
                try:
                    _, _, bits, msg = run_enum(blk, Wk, ck, tau, pout, a.threads, stream=True)
                    dt = np.dtype([('mask', '<u8'), ('key', '<i2', (D,))])
                    n = os.path.getsize(pout) // dt.itemsize
                    mm = np.memmap(pout, dtype=dt, mode='r')
                    masks = np.empty(n, dtype=np.uint64)
                    s1 = np.zeros(D); s2 = np.zeros(D)
                    pick = np.sort(srng.choice(n, min(n, nsamp), replace=False))
                    samp = np.empty((len(pick), D), dtype=np.int16); filled = 0
                    for i in range(0, n, chunk):
                        part = mm[i:i + chunk]; keys = np.ascontiguousarray(part['key'])
                        masks[i:i + len(part)] = part['mask']
                        fk.write(keys.tobytes())
                        X = keys.astype(float); s1 += X.sum(0); s2 += (X * X).sum(0)
                        absmax = np.maximum(absmax, np.abs(keys.astype(np.int64)).max(0))
                        sel = pick[(pick >= i) & (pick < i + len(part))] - i
                        samp[filled:filled + len(sel)] = keys[sel]; filled += len(sel)
                    del mm
                finally:
                    if os.path.exists(pout): os.remove(pout)
                m1 = s1 / max(n, 1); varsum += s2 / max(n, 1) - m1 * m1
                sizes.append(n); samples.append(samp); recs.append(('enum', (masks, bits)))
                blocks[k] = (blk, None, None)
                log(f'block {k}: {msg} (tau {tau:.2f}, rank {rk})')
        with open(path, 'r+b') as fk: fk.write(struct.pack('4i', *sizes))
        est, rank, off = gauss_estimate(samples, t, sizes=sizes)
        log(f'Gaussian estimate of the class: 10^{est0:.2f}; inside the kept lists {sizes}: 10^{est:.2f} (rank {rank}, null offset {off:.2f})')
        fits = absmax <= 1000
        S_idx = [int(x) for x in a.S.split(',')] if a.S else sorted([int(c) for c in np.argsort(-varsum) if fits[c]][:5])
        extra = sum(1 for _ in open(a.skip_sigmas)) if a.skip_sigmas else 0
        sigmas, info = choose_sigmas(samples, S_idx, np.rint(t).astype(np.int64), a.nsig + extra, a.mcpairs, rng, sizes=sizes)
        log(f'split values: {len(sigmas)} chosen from {info["candidates"]} sampled; expected share of the solutions they hold '
            f'{info["share"]:.3g}, predicted pairs AB {info["ab"]:.3g} CD {info["cd"]:.3g}')
        sigmas = skip_done(a, sigmas, info, log)
        del samples
        found = run_kt4(a, S, gens, path, S_idx, sigmas, recs, blocks, Pm0, log)
    finally:
        for pth in (path, path + '.sig'):
            if os.path.exists(pth): os.remove(pth)
    log(f'{len(found)} verified solutions')


def run_kt4(a, S, gens, path, S_idx, sigmas, recs, blocks, Pm0, log):
    """Run kt4 on the list file at path with the given split values; rebuild, check and record every solution."""
    found = []
    env = dict(os.environ)
    if a.threads: env['OMP_NUM_THREADS'] = str(a.threads)
    spath = path + '.sig'
    np.savetxt(spath, sigmas, fmt='%d')
    log(f'kt4: merge classes {S_idx}, {len(sigmas)} split values')
    p = subprocess.Popen([W.bin('kt4'), path, ','.join(map(str, S_idx)), '@' + spath, str(a.K), str(int(a.maxpairs))],
                         stdout=subprocess.PIPE, text=True, env=env)
    for line in p.stdout:
        if line.startswith('TOTAL'): log(f'kt4: {line.strip()}'); continue
        if not line.startswith('SOL'): continue
        sol = list(map(int, line.split()[1:]))
        ids = []
        for k in range(4):
            kind, data = recs[k]
            if kind == 'full': ids.append(Pm0[data[sol[k]]].astype(np.int64))
            else: ids.append(decode(blocks[k][0], data[0][sol[k]], data[1]))
        seqs = [S.sequence(ids[k]) for k in range(4)]
        paf_ok, skew_ok, inv_ok = check_hit(S, seqs, gens)
        ok = paf_ok and skew_ok
        log(f'hit: autocorrelations {"ok" if ok else "WRONG"}, row sums {[int(x.sum()) for x in seqs]}')
        if ok:
            ok = inv_ok
            log(f'hit: invariant under H = <{gens}>: {ok}')
        if ok:
            rec = dict(n=a.v, gens=gens, pattern=[int(x.sum()) for x in seqs], hadamard=True, skew=S.skew,
                       method=f'cp_lowe q={a.q} e={a.e} kgen={a.kgen} twist={a.twist} config={a.config} L={a.L} seed={a.seed}'
                              + (f' field={a.field}' if a.field else ''),
                       seqs=[''.join('+' if z > 0 else '-' for z in x) for x in seqs])
            if a.field: rec.update(group=S.mods, automorphisms=gens)
            with open(W.hits() if S.skew else W.hits('had_hits.jsonl'), 'a') as fh: fh.write(json.dumps(rec) + '\n')
            found.append(rec)
    p.wait()
    return found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('v', type=int); ap.add_argument('--q', type=int, required=True); ap.add_argument('--e', type=int, required=True)
    ap.add_argument('--plain', action='store_true'); ap.add_argument('--config', type=int, default=0)
    ap.add_argument('--kgen', type=int, default=1); ap.add_argument('--twist', type=int, default=0)
    ap.add_argument('--groups', type=int, default=4, help='class groups used as positions when v is prime (q = v)')
    ap.add_argument('--sample-stage1', type=int, default=0, help='sample this many stage-1 quadruples (seeded) instead of listing all')
    ap.add_argument('--L', default='2e7', help='members kept per non-skew block: one value, or one per block 1,2,3 (comma-separated)')
    ap.add_argument('--keepA', type=float, default=0, help='members kept in block 0 (0: all, when it is enumerated in full)')
    ap.add_argument('--samples', type=int, default=200000); ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--S', default='', help='merge classes for kt4 (default: the five of largest spread)')
    ap.add_argument('--nsig', type=int, default=200); ap.add_argument('--maxpairs', type=float, default=3e8)
    ap.add_argument('--mcpairs', type=int, default=4000000, help='random pairs per side for ranking the split values')
    ap.add_argument('--threads', type=int); ap.add_argument('--estimate-only', action='store_true')
    ap.add_argument('--K', type=int, default=10)
    ap.add_argument('--field', default='', help='p,k: the fibre is GF(p^k) (q = p^k; cp_field.FieldStructure)')
    ap.add_argument('--stream', action='store_true', help='stream the lists to the kt4 file (keeps only samples, masks)')
    ap.add_argument('--skip-sigmas', default='', help='file of split values already searched (kt4 log lines or tuples): resume')
    a = ap.parse_args()
    t0 = time.time(); log = lambda s: print(f'{s} [{time.time() - t0:.0f}s]', flush=True)
    Ls = [float(x) for x in a.L.split(',')]
    Lk = lambda k: Ls[0] if len(Ls) == 1 else Ls[k - 1 if len(Ls) == 3 else k]
    if a.field:
        import cp_field as FF
        pf, kf = (int(x) for x in a.field.split(','))
        assert pf ** kf == a.q, '--field p,k must give q'
        S = FF.FieldStructure(a.v, pf, kf, a.e, skew=not a.plain, kgen=a.kgen)
    else:
        S = CP.Structure(a.v, a.q, a.e, skew=not a.plain, kgen=a.kgen)
    O = Orbits(S, a.twist)
    if S.m == 1:
        quads = stage1_prime(S)
        assert quads, 'no stage-1 class'
    elif a.sample_stage1:
        quads = stage1_sample(S, a.sample_stage1, np.random.default_rng(a.seed), O=O)
        assert quads, 'no stage-1 quadruple found by sampling'
    else:
        quads = stage1_cached(S, a.plain)
        if O.tw: quads = [qd for qd in quads if all(admissible(S, O, np.array([y_]))[0] for y_ in qd)]
    y = quads[a.config]
    gens = generators(S, O)
    log(f'v={a.v} = {S.m}*{a.q}, e={a.e}, F={S.F}, K_m={S.K} (twist {O.tw}), H=<{gens}>, {"plain" if a.plain else "skew"}, '
        f'config {a.config} of {len(quads)}: y = {y}')
    rng = np.random.default_rng(a.seed)
    P = prepare_prime(S, y, a.samples, rng, a.groups) if S.m == 1 else prepare(S, O, y, a.samples, rng)
    blocks, full, mus, sigs, Pm0, D = P['blocks'], P['full'], P['mus'], P['sigs'], P['Pm0'], P['D']
    t = np.zeros(D)
    Sig = sum(sigs); Sp = np.linalg.pinv(Sig, rcond=1e-10); mu = sum(mus)
    est0 = P['est0']
    log(f'key dimension {D}, rank {P["rank0"]}; per-block sizes {P["sizes"]}')
    if a.estimate_only:
        est, rank, off = estimate_in_lists(P, Lk)
        log(f'Gaussian estimate of the class: 10^{est0:.2f}; inside low-energy lists of about {a.L}: 10^{est:.2f} (rank {rank}, null offset {off:.2f})')
        return
    if a.stream:
        stream_search(a, S, P, gens, t, Sp, mu, Lk, est0, log, rng)
        return
    lists, recs = [], []
    for k in range(4):
        ck = mus[k] + sigs[k] @ Sp @ (t - mu)
        Wk, rk = whitening(sigs[k])
        if blocks[k] is None:
            K0 = full[k]
            if a.keepA and a.keepA < len(K0):
                en = (((K0 - ck) @ Wk.T) ** 2).sum(1); keep = np.argsort(en)[:int(a.keepA)]
                lists.append(K0[keep]); recs.append(('full', keep))
            else:
                lists.append(K0); recs.append(('full', np.arange(len(K0))))
            log(f'block {k}: {len(lists[-1])} of {len(K0)} members (skew block)')
            continue
        blk, I, K = blocks[k]
        en = (((K - ck) @ Wk.T) ** 2).sum(1)
        frac = min(1.0, Lk(k) / blk.size)
        tau = float(np.quantile(en, frac)) if frac < 1 else 3.0e38
        fd, pout = tempfile.mkstemp(suffix=f'.b{k}', dir=W.tmpdir()); os.close(fd)
        try:
            masks, Kk, bits, msg = run_enum(blk, Wk, ck, tau, pout, a.threads)
        finally:
            os.remove(pout)
        lists.append(Kk); recs.append(('enum', (masks, bits)))
        blocks[k] = (blk, None, None)
        log(f'block {k}: {msg} (tau {tau:.2f}, rank {rk})')
    est, rank, off = gauss_estimate(lists, t)
    log(f'Gaussian estimate of the class: 10^{est0:.2f}; inside the kept lists {[len(l) for l in lists]}: 10^{est:.2f} (rank {rank}, null offset {off:.2f})')
    spread = sum(colvar(l) for l in lists)
    fits = np.max([np.maximum(l.max(0).astype(np.int64), -l.min(0).astype(np.int64)) for l in lists], axis=0) <= 1000   # kt4 packs the merge classes in 12 bits
    S_idx = [int(x) for x in a.S.split(',')] if a.S else sorted([int(c) for c in np.argsort(-spread) if fits[c]][:5])
    extra = sum(1 for _ in open(a.skip_sigmas)) if a.skip_sigmas else 0
    sigmas, info = choose_sigmas(lists, S_idx, np.zeros(D, np.int64), a.nsig + extra, a.mcpairs, rng)
    log(f'split values: {len(sigmas)} chosen from {info["candidates"]} sampled; expected share of the solutions they hold '
        f'{info["share"]:.3g}, predicted pairs AB {info["ab"]:.3g} CD {info["cd"]:.3g}')
    sigmas = skip_done(a, sigmas, info, log)
    fd, path = tempfile.mkstemp(suffix='.kt4', dir=W.tmpdir()); os.close(fd)
    try:
        with open(path, 'wb') as f:
            f.write(struct.pack('5i', *[len(l) for l in lists], D)); f.write(np.zeros(D, np.int16).tobytes())
            for l in lists: f.write(np.ascontiguousarray(l, dtype=np.int16).tobytes())
        del lists                                            # kt4 holds its own copy
        found = run_kt4(a, S, gens, path, S_idx, sigmas, recs, blocks, Pm0, log)
    finally:
        for pth in (path, path + '.sig'):
            if os.path.exists(pth): os.remove(pth)
    log(f'{len(found)} verified solutions')


if __name__ == '__main__':
    main()

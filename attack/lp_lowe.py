"""Legendre pairs on Z_l, l = m q, from H-invariant pairs over Z_m x GF(q), found by a low-energy two-list join.

A Legendre pair (A, B) of odd length l is a pair of +-1 sequences with PAF_A(s) + PAF_B(s) = -2 for every s != 0; then
(sum A)^2 + (sum B)^2 = 2, and after negations both row sums are 1. It gives a Hadamard matrix of order 2l + 2
(Fletcher, Gysin and Seberry, Australas. J. Combin. 23 (2001) 75-86):

        [ -1  -1   1    1  ]
    H = [ -1   1   1   -1  ]          A, B circulant with first rows A, B
        [  1   1   A    B  ]
        [  1  -1   B^T -A^T]

**Structure** (cyclo_product.Structure, plain, with cp_lowe.Orbits): both sequences are invariant under
H = {h = 1 (mod m), h mod q an e-th power} and the multipliers K_m = <kgen> of Z_m, which act as (a, y) -> (k a, g^t y)
(twist t). Over Z_m x GF(q) a sequence is given at each a by f(a) at (a, 0) and its values on the e classes (a, C_c).

**Stage 1** (compressions to Z_m, y(a) = f(a) + F s(a)): K_m-invariant pairs (y_A, y_B) with row sums 1,
PAF_yA(b) + PAF_yB(b) = -2q for b != 0 and sum over the 2m positions of (F s^2 + 2 f s) = 2(m - 1)e (this is
sum y^2 = 2l - 2q + 2).

**Keys** as in cp_lowe: level B (PAF of the compression to GF(q) at the class representatives d_c, target -2m) and
level C (PAF at (b, d_c), 1 <= b <= m/2, target -2). Stage 1 + B + C is equivalent to the Legendre condition.

**Low-energy lists.** With mu_k, Sigma_k the mean and covariance of the keys of block k and Sigma = Sigma_A + Sigma_B,
the key of A in a solution is Gaussian with mean c_A = mu_A + Sigma_A Sigma^-1 (t - mu_A - mu_B) and covariance
Sigma_A - Sigma_A Sigma^-1 Sigma_A, which is Sigma_A / 2 when Sigma_A = Sigma_B; then the whitened offsets z_A, z_B of
the two keys of a solution coincide up to sign and |z|^2 is about rank/2 against rank for a random member. Every member
of both blocks is enumerated by cp_enum.c and kept when |z|^2 <= tau (tau set by the list size L), so a list of a small
fraction of a block still holds a large share of the solutions. The two lists are joined exactly by a linear 64-bit
hash of the keys, and every candidate is rebuilt and checked.

  python3 attack/lp_lowe.py l --q Q --e E [--kgen k] [--twist t] [--config C] [--L 3e7] [--samples 2e5] [--seed 1]
                             [--threads T] [--estimate-only] [--list] [--maxopts 4e6]
      hits -> $HADAMARD_WORK/lp_hits.jsonl
"""
import argparse, itertools, json, math, os, struct, subprocess, sys, tempfile, time
import numpy as np
from math import gcd
from pathlib import Path as _P
sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
import cyclo_product as CP
import cp_lowe as CL

T0 = time.time()


def log(msg):
    print(f'{msg} [{time.time() - T0:.0f}s]', flush=True)


def paf(x):
    x = np.asarray(x, dtype=float)
    return np.rint(np.fft.irfft(np.abs(np.fft.rfft(x)) ** 2, n=len(x))).astype(np.int64)


def is_legendre_pair(A, B):
    p = paf(A) + paf(B)
    return bool(np.all(p[1:] == -2))


# ---- stage 1 ---------------------------------------------------------------------------------------------------
def lp_options(S, O, maxopts=4e6):
    """Compressions y on Z_m of one block of a Legendre pair: constant on the K_m-orbits, entries f + F s with a
    twist-admissible s at every orbit, row sum 1 and spike cost at most 2(m - 1)e. Returns (Yr, C, P, orbits): the
    values at the orbits (one column per orbit), the spike costs and the PAF at b = 1 .. m/2."""
    m, e, F = S.m, S.e, S.F
    budget = 2 * (m - 1) * e
    orbits = [[a for a in range(m) if S.korb[a] == o] for o in range(len(S.kreps))]
    L = [len(o) for o in orbits]
    vals = []                                             # per orbit: (y value, cost) with an admissible pattern
    for orb in orbits:
        r = O.rep_of[orb[0]]
        vs = []
        for s in range(-e, e + 1, 2):
            for f in (1, -1):
                if F * s * s + 2 * f * s > budget: continue
                if any(O.ok(p, r) for p in S.patterns(f, s)): vs.append((f + F * s, len(orb) * (F * s * s + 2 * f * s)))
        vals.append(vs)
    R = len(orbits)
    if any(not vs for vs in vals):                        # an orbit without an admissible pattern
        return np.zeros((0, R), np.int64), np.zeros(0, np.int64), np.zeros((0, m // 2), np.int64), orbits
    reach = [0] * (R + 1)                                 # largest |row sum| the orbits i.. can add
    for i in range(R - 1, -1, -1): reach[i] = reach[i + 1] + L[i] * max(abs(y) for y, _ in vals[i])
    rows, costs = [], []
    cur = [0] * R
    def rec(i, c, rs):
        if abs(1 - rs) > reach[i]: return
        if i == R:
            if rs == 1:
                rows.append(tuple(cur)); costs.append(c)
                if len(rows) > maxopts: raise RuntimeError(f'stage 1: more than {maxopts:.0f} options per block')
            return
        for y, cc in vals[i]:
            if c + cc > budget: continue
            cur[i] = y
            rec(i + 1, c + cc, rs + L[i] * y)
    rec(0, 0, 0)
    Yr = np.array(rows, dtype=np.int64).reshape(-1, R); C = np.array(costs, dtype=np.int64)
    expand = np.zeros(m, dtype=np.int64)
    for i, orb in enumerate(orbits): expand[orb] = i
    P = np.zeros((len(Yr), m // 2), dtype=np.int64)
    for c0 in range(0, len(Yr), 1 << 18):
        Y = Yr[c0:c0 + (1 << 18)][:, expand]
        for b in range(1, m // 2 + 1): P[c0:c0 + len(Y), b - 1] = (Y * np.roll(Y, -b, axis=1)).sum(1)
    return Yr, C, P, orbits


def canon(S, y):
    """Canonical form of a compression under the multipliers u of Z_m (with K_m trivial also the translations)."""
    m = S.m; best = None
    us = [u for u in range(1, m) if gcd(u, m) == 1] if m > 1 else [1]
    for u in us:
        for t in (range(m) if len(S.K) == 1 else [0]):
            z = tuple(int(y[(u * a + t) % m]) for a in range(m))
            if best is None or z < best: best = z
    return best


def lp_stage1(S, O, maxopts=4e6, maxcls=100000):
    """Stage-1 classes: pairs (y_A, y_B) up to the multipliers of Z_m (and translations when K_m is trivial) and the
    order of the two."""
    Yr, C, P, orbits = lp_options(S, O, maxopts)
    m, e, q = S.m, S.e, S.q
    if not len(Yr): return []
    budget = 2 * (m - 1) * e
    rng = np.random.default_rng(7)
    Rv = rng.integers(1, 2 ** 61, size=P.shape[1] + 1)
    with np.errstate(over='ignore'):
        h = (P * Rv[:-1]).sum(1) + C * Rv[-1]
        ht = np.int64((-2 * q * np.ones(P.shape[1], dtype=np.int64) * Rv[:-1]).sum() + budget * Rv[-1])
        need = ht - h
    o = np.argsort(h, kind='stable'); hs = h[o]
    lo = np.searchsorted(hs, need, 'left'); hi = np.searchsorted(hs, need, 'right')
    expand = np.zeros(m, dtype=np.int64)
    for i, orb in enumerate(orbits): expand[orb] = i
    found = {}
    for i in np.flatnonzero(hi > lo):
        for j in o[lo[i]:hi[i]]:
            if j < i or C[i] + C[j] != budget or np.any(P[i] + P[j] != -2 * q): continue
            ya, yb = Yr[i][expand], Yr[j][expand]
            key = tuple(sorted((canon(S, ya), canon(S, yb))))
            if key not in found:
                found[key] = (tuple(int(x) for x in ya), tuple(int(x) for x in yb))
                if len(found) >= maxcls: return [found[k] for k in sorted(found)]
    return [found[k] for k in sorted(found)]


# ---- keys and estimates ----------------------------------------------------------------------------------------
def prepare_lp(S, O, ys, samples, rng):
    pats = set()
    for k in range(2):
        for a in range(S.m):
            r = O.rep_of[a]
            pats |= {O.at(p, a) for p in S.patterns(*CP.fs_of(S, ys[k][r])) if O.ok(p, r)}
    used = sorted(pats)
    idx = np.full(1 << (S.e + 1), -1, dtype=np.int64); idx[used] = np.arange(len(used))
    CC = S.cc_table(used)
    S.bkey = []; S.dkey = S.nb
    blocks, mus, sigs = [], [], []
    for k in range(2):
        blk = CL.Block(S, ys[k], CC, idx, O)
        I = blk.sample(samples, rng); K = blk.keys(I)
        Bv, Cv = CP.level_vectors(S, blk.ids(I[:1000]), idx, CC)
        assert np.array_equal(np.concatenate([Bv, Cv], axis=1), K[:1000]), 'pair tables disagree with level_vectors'
        for j in range(3):                                # cross-check against the PAF of the sequence itself
            x = S.sequence(blk.ids(I[j:j + 1])[0]); p = paf(x)
            shift = lambda b, d: next(u for u in range(S.v) if u % S.m == b % S.m and u % S.q == d % S.q)
            kB = [sum(p[shift(b, S.reps[1 + c])] for b in range(S.m)) for c in range(S.nb)]
            kC = [p[shift(b, S.reps[1 + c])] for b in range(1, S.m // 2 + 1) for c in range(S.e)]
            assert np.array_equal(np.array(kB + kC), K[j]), 'keys disagree with the PAF of the sequence'
        blocks.append((blk, I, K)); mus.append(K.mean(0)); sigs.append(np.cov(K.T.astype(float)))
    D = len(mus[0])
    t = np.array([-2 * S.m] * S.nb + [-2] * (D - S.nb), dtype=float)
    return dict(blocks=blocks, mus=mus, sigs=sigs, D=D, t=t, sizes=[b[0].size for b in blocks], idx=idx, CC=CC)


def gauss2(mus, sigs, t, logN):
    """log10 expected number of pairs with keys summing to t (lattice 4 per coordinate of the range)."""
    mu = mus[0] + mus[1]; Sig = sigs[0] + sigs[1]
    ev, U = np.linalg.eigh(Sig); keep = ev > 1e-9 * ev.max()
    z = U[:, keep].T @ (t - mu)
    off = float(np.abs(U[:, ~keep].T @ (t - mu)).max()) if (~keep).any() else 0.0
    dens = sum(-0.5 * math.log10(2 * math.pi * x) for x in ev[keep]) - 0.5 * float((z * z / ev[keep]).sum()) / math.log(10)
    return logN + dens + int(keep.sum()) * math.log10(4), int(keep.sum()), off


def filters(P):
    """Conditional means c_k, whitenings W_k and the conditional law of the A key given a solution."""
    mus, sigs, t = P['mus'], P['sigs'], P['t']
    Sig = sigs[0] + sigs[1]; Sp = np.linalg.pinv(Sig, rcond=1e-10); mu = mus[0] + mus[1]
    cs = [mus[k] + sigs[k] @ Sp @ (t - mu) for k in range(2)]
    Ws = [CL.whitening(sigs[k])[0] for k in range(2)]
    condcov = sigs[0] - sigs[0] @ Sp @ sigs[0]
    return cs, Ws, condcov


def taus_for(P, L, cs, Ws):
    out = []
    for k in range(2):
        blk, I, K = P['blocks'][k]
        en = (((K - cs[k]) @ Ws[k].T) ** 2).sum(1)
        frac = min(1.0, L[k] / blk.size)
        out.append(float(np.quantile(en, frac)) if frac < 1 else 3.0e38)
    return out


def kept_share(P, cs, Ws, condcov, taus, n=200000, rng=None):
    """Share of the solutions (Gaussian model) whose two keys both pass the energy filters."""
    rng = rng or np.random.default_rng(3)
    ev, U = np.linalg.eigh(condcov); ev = np.clip(ev, 0, None)
    Z = rng.standard_normal((n, len(ev))) * np.sqrt(ev)
    KA = cs[0] + Z @ U.T
    KB = P['t'] - KA
    eA = (((KA - cs[0]) @ Ws[0].T) ** 2).sum(1); eB = (((KB - cs[1]) @ Ws[1].T) ** 2).sum(1)
    return float(np.mean((eA <= taus[0]) & (eB <= taus[1])))


# ---- enumeration and join ----------------------------------------------------------------------------------------
def enum_hashes(blk, Wk, ck, tau, Rv, threads, chunk=1 << 22):
    """Run cp_enum on the block; returns (masks, 64-bit key hashes, bits, message)."""
    fd, pout = tempfile.mkstemp(suffix='.lpb', dir=W.tmpdir()); os.close(fd)
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
        env = dict(os.environ)
        if threads: env['OMP_NUM_THREADS'] = str(threads)
        r = subprocess.run([W.bin('cp_enum'), path_in, pout], capture_output=True, text=True, env=env, check=True)
        os.remove(path_in)
        dt = np.dtype([('mask', '<u8'), ('key', '<i2', (D,))])
        n = os.path.getsize(pout) // dt.itemsize
        masks = np.empty(n, dtype=np.uint64); hs = np.empty(n, dtype=np.int64)
        mm = np.memmap(pout, dtype=dt, mode='r')
        with np.errstate(over='ignore'):
            for i in range(0, n, chunk):
                part = mm[i:i + chunk]
                masks[i:i + len(part)] = part['mask']
                hs[i:i + len(part)] = (part['key'].astype(np.int64) * Rv).sum(1)
        del mm
        return masks, hs, bits, r.stdout.strip()
    finally:
        for p in (path_in, pout):
            if os.path.exists(p): os.remove(p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('l', type=int); ap.add_argument('--q', type=int, required=True); ap.add_argument('--e', type=int, required=True)
    ap.add_argument('--kgen', type=int, default=1); ap.add_argument('--twist', type=int, default=0)
    ap.add_argument('--config', type=int, default=0); ap.add_argument('--L', default='3e7')
    ap.add_argument('--samples', type=int, default=200000); ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--threads', type=int, default=3); ap.add_argument('--maxopts', type=float, default=4e6)
    ap.add_argument('--estimate-only', action='store_true'); ap.add_argument('--list', action='store_true')
    ap.add_argument('--K', type=int, default=3, help='stop after this many verified pairs')
    a = ap.parse_args()
    S = CP.Structure(a.l, a.q, a.e, skew=False, kgen=a.kgen)
    O = CL.Orbits(S, a.twist)
    classes = lp_stage1(S, O, a.maxopts)
    log(f'l={a.l} = {S.m}*{a.q}, e={a.e}, F={S.F}, K_m={S.K} (twist {O.tw}): {len(classes)} stage-1 classes')
    if a.list:
        for i, c in enumerate(classes): print(i, c)
        return
    if not classes: return
    ys = classes[a.config]
    log(f'class {a.config}: y_A = {ys[0]}, y_B = {ys[1]}')
    rng = np.random.default_rng(a.seed)
    P = prepare_lp(S, O, ys, a.samples, rng)
    sizes = P['sizes']; logN = sum(math.log10(x) for x in sizes)
    est, rank, off = gauss2(P['mus'], P['sigs'], P['t'], logN)
    Ls = [float(x) for x in a.L.split(',')]; Ls = Ls * 2 if len(Ls) == 1 else Ls
    cs, Ws, condcov = filters(P)
    taus = taus_for(P, Ls, cs, Ws)
    share = kept_share(P, cs, Ws, condcov, taus)
    log(f'key dimension {P["D"]}, rank {rank} (null offset {off:.2f}); block sizes {sizes}; Gaussian estimate of the class '
        f'10^{est:.2f}; lists of {Ls}: share of the solutions kept {share:.3g} -> 10^{est + math.log10(max(share, 1e-300)):.2f}')
    if a.estimate_only: return
    Rv = np.random.default_rng(11).integers(1, 2 ** 62, size=P['D'])
    tgt = P['t'].astype(np.int64)
    with np.errstate(over='ignore'):
        ht = int((tgt * Rv).sum())
    res = []
    for k in range(2):
        blk = P['blocks'][k][0]
        masks, hs, bits, msg = enum_hashes(blk, Ws[k], cs[k], taus[k], Rv, a.threads)
        res.append((masks, hs, bits)); log(f'block {"AB"[k]}: {msg} (tau {taus[k]:.2f})')
        P['blocks'][k] = (blk, None, None)
    (mA, hA, bA), (mB, hB, bB) = res
    o = np.argsort(hB, kind='stable'); hBs = hB[o]
    with np.errstate(over='ignore'):
        need = np.int64(ht) - hA
    lo = np.searchsorted(hBs, need, 'left'); hi = np.searchsorted(hBs, need, 'right')
    cand = np.flatnonzero(hi > lo)
    log(f'join: {len(cand)} candidate members of A')
    found = 0
    for i in cand:
        for j in o[lo[i]:hi[i]]:
            A = S.sequence(CL.decode(P['blocks'][0][0], mA[i], bA)); B = S.sequence(CL.decode(P['blocks'][1][0], mB[j], bB))
            ok = is_legendre_pair(A, B) and A.sum() == 1 and B.sum() == 1
            log(f'pair: Legendre {"ok" if ok else "WRONG"}')
            if ok:
                rec = dict(l=a.l, order=2 * a.l + 2, q=a.q, e=a.e, kgen=a.kgen, twist=a.twist, config=a.config, L=a.L, seed=a.seed,
                           A=''.join('+' if x > 0 else '-' for x in A), B=''.join('+' if x > 0 else '-' for x in B))
                with open(W.WORK / 'lp_hits.jsonl', 'a') as fh: fh.write(json.dumps(rec) + '\n')
                found += 1
                if found >= a.K: break
        if found >= a.K: break
    log(f'{found} verified Legendre pairs')


if __name__ == '__main__':
    main()

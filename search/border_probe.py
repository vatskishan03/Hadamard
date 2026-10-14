"""Border problems of bordered Goethals-Seidel arrays with Q = Z_q and a character set T (NOTES section 22).

Notation of bgs.py: four +-1 sequences on G with sum_i PAF_i(g) = -Gamma(g + K) (g != 0), Gamma = 4 sum_{chi in T} chi
on Q = G/K, compressions y_i on Q and Mc their Goethals-Seidel array over Q (bordered.gs_array's layout). The border
(X, Rh, Ch) needs X Hadamard of order m = 4t, Rh Rh^T = 4q I, Ch Ch^T = I_4 (x) circ(Gamma) and X Ch^T + Rh Mc^T = 0.

Reformulation used here. The rows of Rh and the columns of Ch are +-1 vectors of V_T (on each of the four blocks,
Fourier support in T), and both form orthogonal bases of V_T. Mc maps V_T to itself with Mc Mc^T = 4t there, so for a
fixed basis Rh the condition reads Ch = -(1/(4t)) Mc Rh^T X: the columns of X are sign vectors eps with
|Mc Rh^T eps| = 4t entrywise (the set S_R), and a border exists for Rh iff S_R holds m mutually orthogonal vectors.
Only the components of the y_i at the characters in T enter.

  python3 search/border_probe.py validate q r e D [--sols N]
      G = Z_q x GF(r) (r small, positions Z_q, cyclotomic classes of index e; K = 0 x GF(r)): exact search for
      sequences with the bordered targets, a border for their compressions, assembly, and the check H H^T = 4n I
      (n = q r + t). D lists the orders of the characters in T (bgs.gamma_vec).
"""
import argparse, itertools, sys
from math import gcd
import numpy as np
from pathlib import Path as _P
sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import spence_lib as L
import bordered as BD
import bgs


def chars_of_orders(q, D):
    """Indices j of the characters chi_j(a) = w^(ja) of Z_q whose order q/gcd(j, q) lies in D."""
    return [j for j in range(q) if q // gcd(j, q) in D]


def block_vectors(q, T):
    """+-1 vectors on Z_q with Fourier support in T."""
    w = np.exp(2j * np.pi / q)
    F = np.array([[w ** (j * a) for a in range(q)] for j in range(q)])
    notT = [j for j in range(q) if j not in T]
    out = []
    for bits in range(1 << q):
        f = np.array([1 - 2 * ((bits >> a) & 1) for a in range(q)])
        if not notT or np.all(np.abs(F[notT] @ f) < 1e-9): out.append(f)
    return np.array(out)


def all_vectors(Vb, nblocks=4):
    idx = np.array(list(itertools.product(range(len(Vb)), repeat=nblocks)))
    return np.concatenate([Vb[idx[:, i]] for i in range(nblocks)], axis=1).astype(np.int8)


def random_basis(ALL, m, rng, tries=500):
    """A random orthogonal basis of V_T made of +-1 vectors: pick, keep the candidates orthogonal to it, repeat."""
    for _ in range(tries):
        cand = np.arange(len(ALL)); chosen = []
        while len(chosen) < m and len(cand):
            v = ALL[cand[rng.integers(len(cand))]].astype(np.int32)
            chosen.append(v)
            cand = cand[(ALL[cand].astype(np.int32) @ v) == 0]
        if len(chosen) == m: return np.array(chosen)
    return None


def S_R(W, t4, m, chunk=1 << 15):
    """Sign vectors eps (eps_0 = +1) with |W eps| = t4 entrywise."""
    out = []
    for c0 in range(0, 1 << (m - 1), chunk):
        idx = np.arange(c0, min(1 << (m - 1), c0 + chunk), dtype=np.int64)
        E = 1 - 2 * ((idx[:, None] >> np.arange(m - 1)[None, :]) & 1)
        E = np.concatenate([np.ones((len(idx), 1), np.int64), E], axis=1)
        ok = np.all(np.abs(E @ W.T) == t4, axis=1)
        if ok.any(): out.append(E[ok])
    return np.concatenate(out) if out else np.zeros((0, m), np.int64)


def ortho_clique(S, m, limit=10 ** 6):
    """m mutually orthogonal rows of S (depth-first search on bitsets), or None."""
    n = len(S)
    if n < m: return None
    G = (S @ S.T) == 0
    adj = [sum(1 << j for j in range(n) if G[i, j]) for i in range(n)]
    best, nodes = [None], [0]
    def rec(cl, cand):
        nodes[0] += 1
        if nodes[0] > limit: return True
        if len(cl) == m: best[0] = list(cl); return True
        if bin(cand).count('1') < m - len(cl): return False
        while cand:
            i = cand.bit_length() - 1; cand &= ~(1 << i)
            if rec(cl + [i], cand & adj[i]): return True
        return False
    rec([], (1 << n) - 1)
    return best[0]


def probe(q, T, ys, nR, rng):
    """Border for compressions ys (four vectors on Z_q) from up to nR random bases: (X, Rh, Ch) or None."""
    t = len(T); m = 4 * t
    ALL = all_vectors(block_vectors(q, T))
    Mc = BD.gs_array(L.AbGroup([q]), [np.asarray(y, float) for y in ys])
    for _ in range(nR):
        R = random_basis(ALL, m, rng)
        if R is None: return None
        W = Mc @ R.T; Wi = np.rint(W).astype(np.int64)
        if not np.allclose(W, Wi): return None
        S = S_R(Wi, 4 * t, m)
        if len(S) >= m:
            cl = ortho_clique(S, m, limit=200000)
            if cl:
                X = S[cl].T
                C = -(Wi @ X) // (4 * t)
                assert np.array_equal(X @ X.T, m * np.eye(m, dtype=np.int64))
                assert np.array_equal(R @ R.T, 4 * q * np.eye(m, dtype=np.int64))
                assert np.array_equal(X @ C.T + R @ np.rint(Mc).astype(np.int64).T, np.zeros((m, 4 * q), np.int64)) or \
                    np.allclose(X @ C.T + R @ Mc.T, 0)
                return X, R, C
    return None


def validate(q, r, e, D, nsol, seed=1):
    """Exhaustive sequences on Z_q x GF(r) (cyclotomic classes of index e), bordered, assembled and checked."""
    S = bgs.Structure(q, r, e, q, list(D), 1)
    T = chars_of_orders(q, D)
    print(S.describe(), 'T =', T, 'Gamma =', S.Gamma.tolist(), flush=True)
    v = q * r
    npat = 1 << (e + 1)
    pats = np.array([[S.pattern(p)[0]] + S.pattern(p)[1] for p in range(npat)])
    tt = np.arange(v); col = np.where(tt % r == 0, 0, 1 + S.cls[tt % r])
    ids = np.array(list(itertools.product(range(npat), repeat=q)))
    X = pats[ids[:, tt % q], col[None, :]]
    P = np.rint(np.fft.ifft(np.abs(np.fft.fft(X.astype(float), axis=1)) ** 2, axis=1).real).astype(np.int64)[:, 1:]
    tgt = np.array([S.target(g) for g in range(1, v)])
    rng = np.random.default_rng(seed)
    Rk = rng.integers(1, 2 ** 61, size=P.shape[1]).astype(np.int64)
    h = (P * Rk).sum(1); n = len(X)
    hs = np.add.outer(h, h).ravel(); order = np.argsort(hs); hsorted = hs[order]
    need = (tgt * Rk).sum() - hs
    pos = np.minimum(np.searchsorted(hsorted, need), len(hsorted) - 1)
    sols = []
    for j in np.flatnonzero(hsorted[pos] == need):
        p = pos[j]
        while p < len(hsorted) and hsorted[p] == need[j]:
            i = order[p]; a, b = divmod(int(i), n); c, d = divmod(int(j), n)
            if np.array_equal(P[a] + P[b] + P[c] + P[d], tgt): sols.append((a, b, c, d))
            p += 1
        if len(sols) >= 20 * nsol: break
    print(f'{n} members per block, {len(sols)} quadruples (capped)', flush=True)
    ok_count = 0
    for a, b, c, d in sols[:nsol]:
        seqs = [X[k] for k in (a, b, c, d)]
        ys = [np.array([s[tt % q == al].sum() for al in range(q)]) for s in seqs]
        res = probe(q, T, ys, 60, rng)
        if res is None:
            print('  compressions', [y.tolist() for y in ys], ': no border found'); continue
        Xb, Rh, Ch = res
        M = BD.gs_array(L.AbGroup([v]), seqs)
        cos = tt % q; mb = Xb.shape[0]
        Rf = np.zeros((mb, 4 * v), np.int64); Cf = np.zeros((4 * v, mb), np.int64)
        for i in range(4):
            Rf[:, i * v:(i + 1) * v] = Rh[:, q * i + cos]; Cf[i * v:(i + 1) * v, :] = Ch[q * i + cos, :]
        Hm = np.block([[Xb, Rf], [Cf, M]])
        ok = np.array_equal(Hm @ Hm.T, Hm.shape[0] * np.eye(Hm.shape[0], dtype=np.int64))
        ok_count += ok
        print(f'  compressions {[y.tolist() for y in ys]}: order {Hm.shape[0]}, H H^T = nI: {ok}', flush=True)
    return ok_count


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    v = sub.add_parser('validate'); v.add_argument('q', type=int); v.add_argument('r', type=int); v.add_argument('e', type=int)
    v.add_argument('D'); v.add_argument('--sols', type=int, default=4); v.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    if a.cmd == 'validate':
        validate(a.q, a.r, a.e, [int(x) for x in a.D.split(',')], a.sols, a.seed)


if __name__ == '__main__':
    main()

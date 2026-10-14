r"""Spence (1975): skew-Hadamard matrices of order 4(1 + q + q^2) of Goethals-Seidel type.

Plane part.  Q = q^2, N = Q^2 + Q + 1 = v w with v = q^2 + q + 1, w = q^2 - q + 1 (coprime).  D = {i mod N : Tr(alpha^i)
= 0} is a Singer (N, Q + 1, 1) difference set (Tr: GF(Q^3) -> GF(Q)); it is fixed by the multiplier q.  Reduce mod v:
f(r) = #{d in D : d = r mod v} lies in {0, 1, 2} (the fibres are orbits of the subgroup of order w, which are arcs), and
    sum_r f(r) f(r + s) = w   for s != 0 mod v            (each difference of D occurs once),
so u = f - 1 has PAF_u(s) = w - 2(Q + 1) + v = 0 off 0.  With D1 = {f = 1}, D2 = {f = 2},
    R = +1 on D1 u D2,  S = +1 on D2   (-1 elsewhere),   R = u + 1_{D1},  S = u - 1_{D1},
    PAF_R + PAF_S = 2 (PAF_u + lambda_{D1}) = 2 lambda_{D1}(s),
which is the constant 2 exactly when D1 is a planar (v, q + 1, 1) difference set (checked below).  A complementary pair
(a skew type, b, PAF_a + PAF_b = -2 off 0, the Szekeres difference sets) over Z_v then gives GS(a, b, R, S), skew of
order 4v.

Usage: spence_plane.py q [--cds FILE]   (FILE: lines 'a +-+...' and 'b +-+...' of length v) [--out cert.json]
"""
import argparse
import os
import json
import random
import sys

import numpy as np
import sympy


def polymulmod(a, b, f, p):
    """a, b: coefficient lists (low first), f monic of degree m"""
    m = len(f) - 1
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] = (r[i + j] + x * y) % p
    for k in range(len(r) - 1, m - 1, -1):
        c = r[k]
        if c:
            for j in range(m + 1):
                r[k - m + j] = (r[k - m + j] - c * f[j]) % p
    r = r[:m] + [0] * max(0, m - len(r))
    return r


def polypowmod(a, n, f, p):
    m = len(f) - 1
    res = [1] + [0] * (m - 1)
    base = a[:] + [0] * (m - len(a))
    while n:
        if n & 1:
            res = polymulmod(res, base, f, p)
        base = polymulmod(base, base, f, p)
        n >>= 1
    return res


def plane_pair(q, seed=1, tries=200):
    fac = sympy.factorint(q)
    (p, e), = fac.items()
    m = 6 * e
    Q = q * q
    N = Q * Q + Q + 1
    v, w = q * q + q + 1, q * q - q + 1
    rng = random.Random(seed)
    for _ in range(tries):
        f = [rng.randrange(p) for _ in range(m)] + [1]
        if f[0] == 0:
            continue
        # companion (multiplication by x) matrix
        M = np.zeros((m, m), dtype=np.int64)
        for j in range(m - 1):
            M[j + 1, j] = 1
        M[:, m - 1] = [(-c) % p for c in f[:m]]
        # Frobenius x -> x^p: column j = x^(p j)
        xp = polypowmod([0, 1], p, f, p)
        F = np.zeros((m, m), dtype=np.int64)
        col = [1] + [0] * (m - 1)
        for j in range(m):
            F[:, j] = col
            col = polymulmod(col, xp, f, p)
        FQ = np.linalg.matrix_power(F, 1) % p
        for _ in range(2 * e - 1):
            FQ = (FQ @ F) % p
        T = (np.eye(m, dtype=np.int64) + FQ + (FQ @ FQ) % p) % p
        # X_i = M^i e0, i < N, by doubling
        X = np.zeros((N, m), dtype=np.int64)
        X[0, 0] = 1
        k, Mk = 1, M.copy()
        while k < N:
            nk = min(k, N - k)
            X[k:k + nk] = (X[:nk] @ Mk.T) % p
            k += nk
            Mk = (Mk @ Mk) % p
        t = (X @ T.T) % p
        D = np.nonzero(~t.any(axis=1))[0]
        if len(D) != Q + 1:
            continue
        diffs = (D[:, None] - D[None, :]) % N
        cnt = np.bincount(diffs.ravel(), minlength=N)
        if cnt[0] != Q + 1 or np.any(cnt[1:] != 1):
            continue
        if set(((q * D) % N).tolist()) != set(D.tolist()):
            raise RuntimeError('D is not fixed by q')
        fr = np.bincount(D % v, minlength=v)
        assert fr.max() <= 2, fr.max()
        D1 = np.nonzero(fr == 1)[0]
        D2 = np.nonzero(fr == 2)[0]
        dd = (D1[:, None] - D1[None, :]) % v
        c1 = np.bincount(dd.ravel(), minlength=v)
        planar = len(D1) == q + 1 and c1[0] == q + 1 and np.all(c1[1:] == 1)
        R = -np.ones(v, dtype=np.int64)
        R[fr >= 1] = 1
        S = -np.ones(v, dtype=np.int64)
        S[fr == 2] = 1
        return dict(q=q, v=v, N=N, f=f, D=D, D1=D1, D2=D2, R=R, S=S, planar=planar)
    raise RuntimeError('no suitable polynomial found')


def paf(x):
    X = np.fft.fft(x)
    return np.rint(np.fft.ifft(X * np.conj(X)).real).astype(np.int64)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('q', type=int)
    ap.add_argument('--cds')
    ap.add_argument('--out')
    a = ap.parse_args()
    P = plane_pair(a.q)
    v = P['v']
    ps = paf(P['R']) + paf(P['S'])
    print(f"q={a.q} v={v} N={P['N']} |D1|={len(P['D1'])} |D2|={len(P['D2'])} planar={P['planar']} "
          f"PAF_R+PAF_S off 0: {sorted(set(ps[1:].tolist()))}")
    if not a.cds:
        return
    seq = {}
    for ln in open(a.cds):
        k, s = ln.split()
        seq[k] = np.array([1 if c == '+' else -1 for c in s], dtype=np.int64)
    A, B = seq['a'], seq['b']
    assert len(A) == v and len(B) == v
    pab = paf(A) + paf(B)
    print('PAF_a + PAF_b off 0:', sorted(set(pab[1:].tolist())), 'a skew type:',
          bool(A[0] == 1 and np.all(A[1:] == -A[1:][::-1])))
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'verify'))
    import verify
    Hm = verify.goethals_seidel(*map(verify.circulant, (A, B, P['R'], P['S'])))
    n = Hm.shape[0]
    had = bool(np.array_equal(Hm @ Hm.T, n * np.eye(n, dtype=np.int64)))
    skew = bool(np.array_equal(Hm + Hm.T, 2 * np.eye(n, dtype=np.int64)))
    print(f'order {n}: Hadamard={had} skew={skew}')
    if a.out and had and skew:
        cert = {'construction': 'goethals-seidel', 'v': v, 'q': a.q,
                'sequences': [''.join('+' if x > 0 else '-' for x in s) for s in (A, B, P['R'], P['S'])]}
        json.dump(cert, open(a.out, 'w'))
        print('wrote', a.out)


if __name__ == '__main__':
    main()

"""octic.py -- exact arithmetic for the octic ring Z[x]/(x^{4v}+1) = Z[w][Z_v]  (v odd, w = primitive 8th root of unity); a module (see verify_all.py).
Z[w] element = integer 4-vector (c0,c1,c2,c3) = c0 + c1 w + c2 w^2 + c3 w^3, w^4 = -1.
Isomorphism  phi: x -> zeta * z,  zeta = w^u, u = v^{-1} mod 8 (= v mod 8), so zeta^v = w.
Inverse: w -> x^v, z -> x^{-8m} with u v = 1 + 8 m.
Sequence a (length 4v, +-1) -> U(g) = sum_{k = g mod v} a_k w^{u k mod 8}.
"""
import numpy as np, itertools, math, random

def zmul(a, b):
    """product in Z[w]; a,b arrays (...,4)"""
    a = np.asarray(a, dtype=np.int64); b = np.asarray(b, dtype=np.int64)
    out = np.zeros(np.broadcast(a, b).shape, dtype=np.int64)
    for i in range(4):
        for j in range(4):
            t = a[..., i] * b[..., j]
            k = i + j
            if k >= 4: out[..., k - 4] -= t
            else: out[..., k] += t
    return out

def zconj(a):
    a = np.asarray(a, dtype=np.int64)
    return np.stack([a[..., 0], -a[..., 3], -a[..., 2], -a[..., 1]], axis=-1)

def zpow_w(k):
    """w^k as 4-vector"""
    k %= 8; s = 1
    if k >= 4: k -= 4; s = -1
    e = np.zeros(4, dtype=np.int64); e[k] = s
    return e

def inv8(v):
    for u in range(1, 8, 2):
        if (u * v) % 8 == 1: return u

def phi(a, v):
    """a: length-4v +-1 (or integer) sequence -> U shape (v,4)"""
    a = np.asarray(a, dtype=np.int64); assert len(a) == 4 * v
    u = inv8(v)
    U = np.zeros((v, 4), dtype=np.int64)
    for k in range(4 * v):
        U[k % v] += a[k] * zpow_w(u * k)
    return U

def gmul(U, V):
    """group ring product over Z_v of Z[w]-valued arrays (v,4)"""
    v = U.shape[0]; out = np.zeros_like(U)
    for g in range(v):
        for h in range(v):
            out[(g + h) % v] += zmul(U[g], V[h])
    return out

def gadj(U):
    """U*(g) = conj(U(-g))"""
    v = U.shape[0]
    return zconj(U[[(-g) % v for g in range(v)]])

def nega_mul(a, b):
    """product of polynomials mod x^n+1 (integer vectors)"""
    n = len(a); out = np.zeros(n, dtype=np.int64)
    for i in range(n):
        for j in range(n):
            k = i + j
            if k >= n: out[k - n] -= a[i] * b[j]
            else: out[k] += a[i] * b[j]
    return out

def nega_adj(a):
    """a(x^{-1}) mod x^n+1: x^{-k} = -x^{n-k}"""
    n = len(a); out = np.zeros(n, dtype=np.int64); out[0] = a[0]
    for k in range(1, n): out[n - k] = -a[k]
    return out

def nega_paf(a, s):
    n = len(a); return int(sum((a[i] * a[i + s] if i + s < n else -a[i] * a[i + s - n]) for i in range(n)))

def alphabet():
    return [np.array([c0, c1, c2, c3]) for c0 in (1, -1) for c1 in (1, -1) for c2 in (1, -1) for c3 in (1, -1)]

def zcomplex(a):
    w = np.exp(1j * np.pi / 4)
    a = np.asarray(a)
    return sum(a[..., k] * w ** k for k in range(4))

def zeta_eval(a, v):
    """R = a(zeta) in Z[w], zeta = w^u  (the row sum / DC value)"""
    u = inv8(v); R = np.zeros(4, dtype=np.int64)
    for k, ak in enumerate(a): R += ak * zpow_w(u * k)
    return R

def psi(U, v):
    """inverse of phi: U (v,4) -> integer sequence of length 4v"""
    n = 4 * v; u = inv8(v); m = (u * v - 1) // 8
    out = np.zeros(n, dtype=np.int64)
    for g in range(v):
        for j in range(4):
            if U[g, j]:
                k = (v * j - 8 * m * g) % (8 * v); s = 1
                if k >= 4 * v: k -= 4 * v; s = -1
                out[k] += s * U[g, j]
    return out

def code_to_c(e):
    return np.array([-1 if (e >> j) & 1 else 1 for j in range(4)], dtype=np.int64)

def expand_codes(codes, orb):
    """codes: string of hex digits (one per orbit) -> U array (v,4)"""
    v = len(orb)
    return np.array([code_to_c(int(codes[orb[g]], 16)) for g in range(v)], dtype=np.int64)

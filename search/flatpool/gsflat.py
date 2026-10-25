#!/usr/bin/env python3
"""Structure tables for flat-pool Goethals-Seidel searches over Z_v, v prime, U <= Z_v^* of odd order.

Orbits: {0} and the m = (v-1)/|U| cosets rho^j U (rho = primitive root; j in Z_m).  Multiplication by rho rotates j.
-1 = rho^(m/2) (m even) rotates by m/2.  C = m/2 classes of d != 0 under <U,-1>, class of rho^e = e mod C.
Spectrum of the sequence (f0; s_0..s_{m-1}):  fhat_c = f0 + sum_j s_j eta_{(j+c) mod m},  c in Z_C,
eta_k = sum_{u in U} exp(2 pi i rho^k u / v);  PSD_c = |fhat_c|^2 (c = class of the frequency rho^c).
PAF(rho^e) = (r^2 + sum_c PSD_c * 2 Re eta_{c+e}) / v,  r = f0 + |U| sum_j s_j (row sum).
Key (integer): kappa_e = (PAF(rho^e) - (v mod 4)) / 4; the four keys must sum to -(v mod 4) in every class.
"""
import math, sys, itertools
import numpy as np
from sympy import primitive_root, factorint


def build(v, ulen, gen=None):
    rho = primitive_root(v)
    assert (v - 1) % ulen == 0 and ulen % 2 == 1
    m = (v - 1) // ulen
    assert m % 2 == 0
    U = sorted({pow(rho, m * t, v) for t in range(ulen)})
    if gen is not None:
        S = {1}; x = 1
        while True:
            x = x * gen % v
            if x == 1: break
            S.add(x)
        assert sorted(S) == U, "generator does not generate the order-%d subgroup" % ulen
    assert (v - 1) not in U
    dlog = np.zeros(v, dtype=np.int64)
    x = 1
    for k in range(v - 1):
        dlog[x] = k; x = x * rho % v
    C = m // 2
    eta = np.zeros(m, dtype=complex)
    for k in range(m):
        rk = pow(rho, k, v)
        eta[k] = sum(np.exp(2j * np.pi * ((rk * u) % v) / v) for u in U)
    return dict(v=v, ulen=ulen, rho=rho, m=m, C=C, U=U, dlog=dlog, eta=eta, pm=v % 4)


def full_seq(S, f0, s):
    v, m, dlog = S['v'], S['m'], S['dlog']
    f = np.zeros(v, dtype=np.int64)
    f[0] = f0
    f[1:] = np.array(s)[dlog[1:] % m]
    return f


def paf_exact(f):
    v = len(f)
    F = np.fft.fft(f.astype(float))
    p = np.fft.ifft(np.abs(F) ** 2).real
    return np.rint(p).astype(np.int64)


def paf_classes_exact(S, f):
    p = paf_exact(f)
    rho, v, C = S['rho'], S['v'], S['C']
    return np.array([p[pow(rho, e, v)] for e in range(C)])


def spectrum(S, f0, s):
    m, C, eta = S['m'], S['C'], S['eta']
    s = np.asarray(s, float)
    fh = np.array([f0 + sum(s[j] * eta[(j + c) % m] for j in range(m)) for c in range(C)])
    return fh


def paf_from_spectrum(S, f0, s):
    m, C, eta, v = S['m'], S['C'], S['eta'], S['v']
    r = f0 + S['ulen'] * int(sum(s))
    psd = np.abs(spectrum(S, f0, s)) ** 2
    th = 2 * eta.real
    return np.array([(r * r + sum(psd[c] * th[(c + e) % m] for c in range(C))) / v for e in range(C)]), psd


def rowsum_options(v, ulen):
    """Pairs (r, f0, S, npos_cosets) with r = f0 + ulen*S, r>0, attainable.  S even (m even)."""
    m = (v - 1) // ulen
    out = {}
    for f0 in (1, -1):
        for Ssum in range(-m, m + 1, 2):
            r = f0 + ulen * Ssum
            if r > 0:
                out.setdefault(r, []).append((f0, Ssum, (m + Ssum) // 2))
    return out


def patterns(v, ulen):
    opts = rowsum_options(v, ulen)
    rs = sorted(opts)
    out = []
    for a, b, c, d in itertools.combinations_with_replacement(rs, 4):
        if a * a + b * b + c * c + d * d == 4 * v:
            out.append((a, b, c, d))
    return out


if __name__ == '__main__':
    v, ulen = int(sys.argv[1]), int(sys.argv[2])
    gen = int(sys.argv[3]) if len(sys.argv) > 3 else None
    S = build(v, ulen, gen)
    print("v", v, "|U|", ulen, "rho", S['rho'], "m", S['m'], "C", S['C'], "patterns", patterns(v, ulen))
    rng = np.random.default_rng(1)
    opts = rowsum_options(v, ulen)
    # check spectrum / PAF formulas against exact computation on random sequences
    for trial in range(3):
        s = rng.choice([-1, 1], size=S['m'])
        f0 = int(rng.choice([-1, 1]))
        f = full_seq(S, f0, s)
        pe = paf_classes_exact(S, f)
        pf, psd = paf_from_spectrum(S, f0, s)
        fh = np.abs(np.fft.fft(f.astype(float))) ** 2
        rho = S['rho']
        psd_ex = np.array([fh[pow(rho, c, v)] for c in range(S['C'])])
        print("maxdiff PAF", np.abs(pe - pf).max(), "PSD", np.abs(psd - psd_ex).max(), "PAF mod 4 ->", set((pe % 4).tolist()))

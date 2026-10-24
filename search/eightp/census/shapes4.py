"""shapes4.py -- shape enumeration and exact filters (identity / row-sum DP in orbit units) for the real plug-in modes a, b, c."""
import sys, json, math, itertools, time
import numpy as np
from math import gcd
from scipy.signal import fftconvolve
import os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'gs'))
from k4 import *


def prime_power(q):
    for p in range(2, q + 1):
        if q % p == 0:
            e = 0; t = q
            while t % p == 0: t //= p; e += 1
            return (p, e) if t == 1 else None
    return None


def odd_prime_powers_dividing(n, rmin=5):
    out = []
    d = 2; f = {}
    t = n
    while d * d <= t:
        while t % d == 0: f[d] = f.get(d, 0) + 1; t //= d
        d += 1
    if t > 1: f[t] = f.get(t, 0) + 1
    for p, a in f.items():
        if p % 2 == 1:
            for e in range(1, a + 1):
                if p ** e >= rmin: out.append(p ** e)
    return sorted(out)


def orbit_data(sh):
    """list of (weight, forced) of the CD orbits"""
    K = sh.kCD
    return [(int(K.w[i]), bool(K.forced[i])) for i in range(K.nu)]


def block_table(orbs, r, Cmax, g, Rmax):
    """bool table A[c/g][rho + Rmax]: achievable (cw, rho) of one block"""
    nC = Cmax // g + 1
    A = np.zeros((nC, 2 * Rmax + 1), dtype=bool)
    A[0, Rmax] = True
    for (w, f) in orbs:
        B = np.zeros_like(A)
        variants = []
        if not f:
            variants += [(0, w), (0, -w)]
        variants += [(((r + 1) * w) // g, r * w), (((r + 1) * w) // g, -r * w), (((r - 3) * w) // g, (r - 2) * w), (((r - 3) * w) // g, -(r - 2) * w)]
        for dc, dr in variants:
            if dc >= nC or abs(dr) > 2 * Rmax: continue
            src = A[:nC - dc]
            if dr >= 0:
                B[dc:, dr:] |= src[:, :2 * Rmax + 1 - dr]
            else:
                B[dc:, :2 * Rmax + 1 + dr] |= src[:, -dr:]
        A = B
    return A


def dp_filter(sh):
    """EXACT (in orbit units): exists an assignment of (nondefect / defect p / defect m, sign) to the orbits of the four blocks with the defect identity
    (r+1) e_p + (r-3) e_m = 4|K| + phi(0) and sum_i rho_i^2 = rowsq (rho_i = sum_k u_i(k) the block sums).  Returns (ok, info)."""
    r = sh.q; N = sh.N; D = sh.Did; R2 = int(sh.rowsq)
    if D < 0: return False, 'identity: negative target'
    orbs = orbit_data(sh)
    if r == 3: pass
    incs = [(r + 1) * w for w, f in orbs] + ([(r - 3) * w for w, f in orbs] if r > 3 else [])
    g = D
    for x in incs:
        if x: g = gcd(g, x)
    if g == 0: g = 1
    if D % g: return False, 'identity: gcd'
    if sum(w for w, f in orbs if f) * (r - 3 if r > 3 else 1) * 4 > D + 4 * N * (r == 3): pass
    rho_max = int(math.isqrt(max(R2, 0)))
    Rmax = min(N + (D // max(r - 3, 1) + 1) * r, 4 * N * r) + 1
    Rmax = min(Rmax, 6000)
    A = block_table(orbs, r, D, g, Rmax)
    if not A.any(): return False, 'identity(forced defects too heavy)'
    nC = D // g + 1
    a1 = A.any(axis=1).astype(float)
    p2 = np.convolve(a1, a1)[:nC] > 0.5
    p4 = np.convolve(p2.astype(float), p2.astype(float)) > 0.5
    if len(p4) <= nC - 1 or not p4[nC - 1]: return False, 'identity(no orbit-unit defect allocation)'
    A2 = np.zeros((nC, R2 + 1), dtype=np.float32)
    rhos = np.arange(-Rmax, Rmax + 1)
    for j in np.flatnonzero(np.abs(rhos) <= rho_max):
        v = int(rhos[j]) ** 2
        if v <= R2: A2[:, v] = np.maximum(A2[:, v], A[:, j])
    if not A2.any(): return False, 'rowsum: no block with |rho| small'
    P = fftconvolve(A2, A2)[:nC, :R2 + 1] > 0.5
    # P[c][v]: pairs; need c1 + c2 = nC-1, v1 + v2 = R2
    Pf = P.astype(np.float32)
    ok = False
    for c in range(nC):
        if P[c].any() and P[nC - 1 - c].any():
            if (P[c] & P[nC - 1 - c][::-1]).any(): ok = True; break
    if ok: return True, 'ok'
    # distinguish: identity alone vs identity + rowsum
    return False, 'rowsum(identity ok, row sums impossible)'


def abelian_orderings(Kmods, qkind):
    """for mode c: reorder the prime-power moduli so that the quotient coordinates come first; returns list of (Kmods_reordered, qspec)"""
    out = []
    if qkind == 'z3':
        idx = [i for i, m in enumerate(Kmods) if m % 3 == 0]
        seen = set()
        for i in idx:
            if Kmods[i] in seen: continue
            seen.add(Kmods[i])
            K2 = [Kmods[i]] + [m for j, m in enumerate(Kmods) if j != i]
            out.append((K2, ([3], [(0,), (1,), (2,)])))
    elif qkind == 'z2z2':
        ev = [i for i, m in enumerate(Kmods) if m % 2 == 0]
        seen = set()
        for a, b in itertools.combinations(ev, 2):
            key = tuple(sorted((Kmods[a], Kmods[b])))
            if key in seen: continue
            seen.add(key)
            K2 = [Kmods[a], Kmods[b]] + [m for j, m in enumerate(Kmods) if j not in (a, b)]
            out.append((K2, ([2, 2], [(1, 0), (0, 1), (1, 1)])))
    elif qkind.startswith('z') and qkind[1:].isdigit():
        q = int(qkind[1:])
        idx = [i for i, m in enumerate(Kmods) if m % q == 0]
        seen = set()
        for i in idx:
            if Kmods[i] in seen: continue
            seen.add(Kmods[i])
            K2 = [Kmods[i]] + [m for j, m in enumerate(Kmods) if j != i]
            out.append((K2, ([q], [(a,) for a in range(q)])))
    return out


def enumerate_shapes(n, modes='abc', Orange=(4, 28), rmax=10 ** 9, qkinds=('z3', 'z2z2'), twist=True, log=None, rmin=5):
    """yield descriptors dict(mode, n, r, K, hg, qspec, qkind, t, orb, forced, ...) for all (r, K, H) with the orbit count in the window"""
    for mode in modes:
        if mode == 'a': tlist = [(0, None)]
        elif mode == 'b': tlist = [(1, 'b')]
        else: tlist = [(3, qk) for qk in qkinds if qk in ('z3', 'z2z2')] + [(int(qk[1:]), qk) for qk in qkinds if qk.startswith('z') and qk[1:].isdigit() and int(qk[1:]) not in (3,) and qk != 'z2z2']
        for (t, qk) in tlist:
            Gsz = n - t
            if Gsz < 6: continue
            if mode == 'c' and qk == 'z2z2' and n % 4 != 3: continue
            for r in odd_prime_powers_dividing(Gsz, rmin):
                if r > rmax: continue
                N = Gsz // r
                if N < 4: continue
                for Km0 in abelian_groups_pp(N):
                    if mode == 'a': variants = [(Km0, None)]
                    elif mode == 'b': variants = [(Km0, ([], [()]))]
                    else: variants = abelian_orderings(Km0, qk)
                    for (Km, qs) in variants:
                        for hg, hs in cyclic_unit_subgroups(Km, twist=twist):
                            no, nf = orbit_count(Km, hg)
                            if no < Orange[0] or no > Orange[1]: continue
                            yield dict(mode=mode, n=n, r=r, K=Km, hg=[[list(u), ph] for u, ph in hg], qspec=qs, qkind=qk, t=t, orb=no, forced=nf, hsz=hs, N=N)

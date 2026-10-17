#!/usr/bin/env python3
"""Admissibility scan for T-matrices of prime order t with multiplier symmetry, for Cooper-Wallis
routes to open orders 4n with n = t*w, w a Williamson order.

T-matrices: four circulant 0/+-1 matrices on Z_t, supports partition Z_t, sum PAF = t*delta.
Row sums: sum r_i^2 = t.  With H <= Z_t^* of order h, supports are unions of H-orbits:
 plain: block containing 0 has r = +-1 + h*j, the others r = h*j.
 T3 (K = <H, mu>, [K:H]=3): (T1, T2, mu T2, mu^2 T2), T1 K-invariant contains 0:
        r1 = +-1 + 3h*j, r = h*j',  r1^2 + 3 r^2 = t.
Classes: <K,-1>-orbits on Z_t minus 0."""
import os
import re, math
from sympy import isprime, divisors

WILL = set(range(1, 34, 2)) | {37, 39, 41, 43, 45, 49, 51, 55, 57, 61, 63}

def open_ns():
    ns = []
    for line in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'FRONTIER.md')):
        m = re.match(r'\|\s*(\d+)\s*\|\s*(\d+)\s*\|', line)
        if m: ns.append(int(m.group(2)))
    return ns

def classes(t, d):
    return (t - 1) // d if d % 2 == 0 else (t - 1) // (2 * d)

def plain(t, h):
    k = (t - 1) // h
    for r0 in range(-(h * k) - 1, h * k + 2):
        if (r0 - 1) % h and (r0 + 1) % h: continue
        rest = t - r0 * r0
        if rest < 0 or rest % (h * h): continue
        m = rest // (h * h)
        # m must be a sum of three squares
        for a in range(int(math.isqrt(m)) + 1):
            for b in range(a, int(math.isqrt(m - a * a)) + 1):
                c2 = m - a * a - b * b
                c = math.isqrt(c2)
                if c * c == c2 and c >= b:
                    return (abs(r0), h * a, h * b, h * c)
    return None

def twist3(t, h):
    d = 3 * h
    for r in range(1, int(math.isqrt(t // 3)) + 1):
        if r % h: continue
        r1s = t - 3 * r * r
        r1 = math.isqrt(r1s)
        if r1 * r1 == r1s and ((r1 - 1) % d == 0 or (r1 + 1) % d == 0):
            return (r1, r)
    return None

rows = []
for n in open_ns():
    for w in sorted(WILL):
        if n % w: continue
        t = n // w
        if not isprime(t) or t < 70: continue
        for h in divisors(t - 1):
            if h == 1 or h == t - 1: continue
            pat = plain(t, h)
            if pat: rows.append((classes(t, h), n, w, t, 'plain', h, h, pat, (t - 1) // h))
            if (t - 1) % (3 * h) == 0:
                pat = twist3(t, h)
                if pat: rows.append((classes(t, 3 * h), n, w, t, 'T3', h, 3 * h, pat, (t - 1) // (3 * h)))
rows.sort()
for r in rows:
    if r[0] <= 40:
        print(f"C={r[0]:3d}  4n={4*r[1]:6d} n={r[1]:5d} = {r[2]:2d} x t={r[3]:5d}  {r[4]:5s} h={r[5]:4d} |K|={r[6]:4d} K-orbits={r[8]:4d} rows={r[7]}")

#!/usr/bin/env python3
r"""Admissibility scan of twisted multiplier structures over Z_p (p prime).

K <= Z_p^* of order d acts on the four Goethals-Seidel blocks through K/H:
  plain : all blocks K-invariant                         sum r_i^2 = 4p
  T3    : [K:H]=3, (A, B, B.mu, B.mu^2), A K-inv, B H-inv  r_A^2 + 3 r_B^2 = 4p
  T2a   : [K:H]=2, (A, A.mu, C, C.mu), A,C H-inv           2 r_A^2 + 2 r_C^2 = 4p
  T2b   : [K:H]=2, (A, A.mu, C, D), A H-inv, C,D K-inv     2 r_A^2 + r_C^2 + r_D^2 = 4p
Row sum of a sequence constant on orbits of size m (k = (p-1)/m orbits) is
  +-1 + m*j,  j = k mod 2, |j| <= k.
Classes: <K,-1>-orbits on Z_p \ {0}.
"""
import sys, math, re
from sympy import isprime, divisors

def rowsums(p, m):
    k = (p - 1) // m
    return sorted({e + m * j for e in (1, -1) for j in range(-k, k + 1, 2)})

def admissible(p, d, kind):
    four = 4 * p
    if kind == 'plain':
        R = [r for r in rowsums(p, d) if r > 0 and r * r <= four]
        sq = {r * r for r in R}
        for a in R:
            for b in R:
                for c in R:
                    rest = four - a*a - b*b - c*c
                    if rest > 0 and rest in sq:
                        return (a, b, c, int(math.isqrt(rest)))
        return None
    if kind == 'T3':
        if d % 3: return None
        RA = [r for r in rowsums(p, d) if r > 0]
        RB = [r for r in rowsums(p, d // 3) if r > 0]
        for a in RA:
            for b in RB:
                if a*a + 3*b*b == four: return (a, b)
        return None
    if d % 2: return None
    h = d // 2
    if kind == 'T2a':
        R = [r for r in rowsums(p, h) if r > 0]
        for a in R:
            for c in R:
                if 2*a*a + 2*c*c == four: return (a, c)
        return None
    if kind == 'T2b':
        RA = [r for r in rowsums(p, h) if r > 0]
        RC = [r for r in rowsums(p, d) if r > 0]
        sq = set(x*x for x in RC)
        for a in RA:
            for c in RC:
                rest = four - 2*a*a - c*c
                if rest > 0 and rest in sq: return (a, c, int(math.isqrt(rest)))
        return None

def classes(p, d):
    minus1_in_K = ((p - 1) // d) % 2 == 0 if False else pow(p - 1, d, p) == 1  # -1 in K iff (-1)^... K = d-th roots
    # K = {x : x^d = 1}; -1 in K iff (-1)^d = 1 iff d even
    minus1_in_K = d % 2 == 0
    return (p - 1) // d if minus1_in_K else (p - 1) // (2 * d)

def open_ns():
    ns = []
    for line in open('../FRONTIER.md'):
        m = re.match(r'\|\s*(\d+)\s*\|\s*(\d+)\s*\|', line)
        if m: ns.append(int(m.group(2)))
    return ns

if __name__ == '__main__':
    maxC = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    rows = []
    for p in open_ns():
        if not isprime(p): continue
        for d in divisors(p - 1):
            if d == 1: continue
            C = classes(p, d)
            if C > maxC: continue
            for kind in ('plain', 'T3', 'T2a', 'T2b'):
                if kind == 'T2a' and d % 2 == 0 and (d // 2) % 2 == 1 and d % 2 == 0:
                    pass
                pat = admissible(p, d, kind)
                if pat:
                    rows.append((C, p, d, kind, pat))
    rows.sort()
    for r in rows:
        print(f"C={r[0]:3d}  p={r[1]:5d}  |K|={r[2]:4d}  {r[3]:5s}  rows={r[4]}")

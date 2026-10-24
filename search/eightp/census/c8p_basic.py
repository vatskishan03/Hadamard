"""c8p_basic.py -- exact census of the plain multiplier structures for the six open orders 8p
(p = 599, 719, 739, 839, 1031, 1087), and the row-sum (DC) filters.  All outputs EXACT.

 (a) real Goethals-Seidel over Z_2p = Z_2 x Z_p, multiplier group H <= Z_p^* of odd order h (h | (p-1)/2):
     orbits = 2 + 2(p-1)/h,  row sums r_X = s0 + s1, s_a = +-1 + 2 h k  (each Z_p half)
     exact conditions from the compression to Z_2:  sum_X (s0^2 + s1^2) = 8p,  sum_X s0 s1 = 0
 (b) negacyclic GS of order 8p = quaternary sequences on Z_p, sum u u* = 4p delta, H of odd order h:
     orbits = 1 + (p-1)/h,  block sum R = u0 + h S, S in Z[i], S = a+bi, |a|+|b| <= c, a+b = c mod 2, c = (p-1)/h
     exact DC condition sum |R_X|^2 = 4p.
"""
import sys, math, itertools
import numpy as np
from sympy import factorint, isprime, divisors

PRIMES = [599, 719, 739, 839, 1031, 1087]


def dc_real_z2p(p, h):
    """exact feasibility of the Z_2-compression + DC conditions of real GS over Z_2 x Z_p with H of odd order h"""
    c = (p - 1) // h            # number of H-orbits on Z_p^*  (even)
    assert c % 2 == 0
    # s = e0 + h*sum_{c orbits} e_j : e0 = +-1, sum in {-c, -c+2, ..., c}
    vals = set()
    for e0 in (1, -1):
        for m in range(-c, c + 1, 2):
            vals.add(e0 + h * m)
    vals = sorted(vals)
    bound = 8 * p
    pairs = {}
    for s0 in vals:
        if s0 * s0 > bound: continue
        for s1 in vals:
            q = s0 * s0 + s1 * s1
            if q > bound: continue
            pairs.setdefault((q, s0 * s1), []).append((s0, s1))
    keys = list(pairs)
    # four pairs summing to (8p, 0)
    two = {}
    for a in keys:
        for b in keys:
            k = (a[0] + b[0], a[1] + b[1])
            if k[0] <= bound: two.setdefault(k, (a, b))
    for k, (a, b) in two.items():
        need = (bound - k[0], -k[1])
        if need in two:
            return True, (a, b, two[need])
    return False, None


def dc_quat_zp(p, h):
    c = (p - 1) // h
    vals = {}
    bound = 4 * p
    u0s = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    # S = a + b i, |a|+|b| <= c, a+b = c mod 2 ; only |R|^2 <= bound matters
    amax = int(math.isqrt(bound)) // h + 2
    for (ur, ui) in u0s:
        for a in range(-min(c, amax), min(c, amax) + 1):
            for b in range(-min(c, amax), min(c, amax) + 1):
                if abs(a) + abs(b) > c or (a + b - c) % 2: continue
                R2 = (ur + h * a) ** 2 + (ui + h * b) ** 2
                if R2 <= bound: vals.setdefault(R2, (ur, ui, a, b))
    ks = sorted(vals)
    two = {}
    for a in ks:
        for b in ks:
            if a + b <= bound: two.setdefault(a + b, (a, b))
    for s, (a, b) in two.items():
        if bound - s in two:
            return True, (a, b, two[bound - s])
    return False, None


def main():
    for p in PRIMES:
        fp = factorint(p - 1)
        print(f"p = {p}  p mod 8 = {p % 8}  p-1 = {dict(fp)}  (p-1)/2 = {(p-1)//2} = {dict(factorint((p-1)//2))}")
        print(f"   2p-1 = {2*p-1} = {dict(factorint(2*p-1))};  2p-2 = {dict(factorint(2*p-2))};  2p-3 = {dict(factorint(2*p-3))};  2p+1 = {dict(factorint(2*p+1))}")
        for h in divisors((p - 1) // 2):
            c = (p - 1) // h
            ok_r, w_r = dc_real_z2p(p, h)
            ok_q, w_q = dc_quat_zp(p, h)
            print(f"   h = {h:4d}: real Z_2p orbits {2 + 2*c:5d} (2^{2+2*c} per block) DC {'live' if ok_r else 'DEAD'};"
                  f"  quaternary Z_p orbits {1 + c:4d} (4^{1+c} = 2^{2+2*c}) DC {'live' if ok_q else 'DEAD'}"
                  + (f"  e.g. |R|^2 {w_q}" if ok_q and h > 1 else ""))
    # largest feasible h for the plain bound
    print()
    print("(h = 1 is the trivial group: 2^{2p} per block, never searchable.)")


if __name__ == '__main__':
    main()

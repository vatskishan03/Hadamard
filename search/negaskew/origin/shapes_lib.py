import itertools
from math import gcd
from sympy import isprime
def subgroups_odd(m):
    U = [u for u in range(1, m) if gcd(u, m) == 1]
    def order(u):
        x, k = u, 1
        while x != 1: x = x * u % m; k += 1
        return k
    odd = [u for u in U if order(u) % 2 == 1 and u != 1]
    def closure(gens):
        H = {1}; fr = [1]
        while fr:
            x = fr.pop()
            for g in gens:
                y = x * g % m
                if y not in H: H.add(y); fr.append(y)
        return frozenset(H)
    subs = {frozenset([1]): [1]}
    for u in odd:
        H = closure([u]); subs.setdefault(H, [u])
    base = list(subs.items())
    for (H1, g1), (H2, g2) in itertools.combinations(base, 2):
        if H1 <= H2 or H2 <= H1: continue
        H = closure(g1 + g2)
        if len(H) % 2 == 0: continue
        subs.setdefault(H, g1 + g2)
    return subs
def orbits(m, H):
    seen = [False] * m; orbs = []
    for k in range(m):
        if not seen[k]:
            o = sorted({k * h % m for h in H})
            for x in o: seen[x] = True
            orbs.append(o)
    return orbs
def primes34(lo, hi):
    return [q for q in range(lo, hi + 1) if isprime(q) and q % 4 == 3]

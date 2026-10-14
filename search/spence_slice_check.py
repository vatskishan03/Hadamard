"""Slice-character admissibility of Spence structures over G = Z_2 x Z_m (m = (n-1)/2 odd).

For H <= Z_m^* acting on Z_m (trivially on Z_2), sequences constant on H-orbits have slice sums
a_i (on {0} x Z_m) and b_i (on {1} x Z_m), each a signed sum of the K-orbit sizes.  The order-2
character of G forces sum_i (a_i - b_i)^2 = 4n, with a_i + b_i = r_i = (2,0,0,0); equivalently
(a_1 - 1)^2 + a_2^2 + a_3^2 + a_4^2 = n.  Returns whether any H (orbits <= MAXO) admits this."""
import math, re, sys
from math import gcd

def subgroups(m):
    units = [u for u in range(1, m) if gcd(u, m) == 1]
    def closure(gens):
        S, todo = {1}, [1]
        while todo:
            x = todo.pop()
            for g in gens:
                y = x * g % m
                if y not in S: S.add(y); todo.append(y)
        return frozenset(S)
    cyc = {}
    for u in units: cyc.setdefault(closure([u]), u)
    subs = set(cyc)
    items = list(cyc.items())
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            subs.add(closure([items[i][1], items[j][1]]))
    return subs

def orbit_sizes(m, H):
    seen, sizes = [False] * m, []
    for t in range(m):
        if not seen[t]:
            O = {t * h % m for h in H}
            for x in O: seen[x] = True
            sizes.append(len(O))
    return sizes

def achievable(sizes):
    S = {0}
    for s in sizes:
        S = {x + s for x in S} | {x - s for x in S}
    return S

def admissible(n, sizes):
    S = achievable(sizes)
    A1 = [a for a in S if (2 - a) in S]
    sq = {}
    for a in S:
        sq.setdefault(a * a, a)
    squares = sorted(sq)
    for a1 in A1:
        rest = n - (a1 - 1) ** 2
        if rest < 0: continue
        for x in squares:
            if x > rest: break
            for y in squares:
                if x + y > rest: break
                z = rest - x - y
                if z in sq: return True, (a1, sq[x], sq[y], sq[z])
    return False, None

if __name__ == "__main__":
    fr = open(sys.argv[1]).read()
    MAXO = int(sys.argv[2]) if len(sys.argv) > 2 else 400
    rows = re.findall(r'^\|\s*(\d+)\s*\|\s*(\d+)\s*\|(.*)$', fr, re.M)
    bad = good = 0
    for N4, n, rest in rows:
        n = int(n)
        cols = [c.strip() for c in rest.split('|')]
        spence_group, spence_w = cols[5], cols[6]
        if spence_group in ('–', '') or n % 4 != 3: continue
        m = (n - 1) // 2
        found = None
        for H in sorted(subgroups(m), key=len, reverse=True):
            sizes = orbit_sizes(m, H)
            if 2 * len(sizes) > MAXO: continue
            ok, wit = admissible(n, sizes)
            if ok: found = (len(H), 2 * len(sizes), wit); break
        if found: good += 1
        else: bad += 1
        print(f"{N4:>6} n={n:<5} FRONTIER Spence {spence_group:22s} logW {spence_w:>5}:  "
              + (f"admissible with |H|={found[0]}, {found[1]} orbits, slice sums {found[2]}" if found
                 else f"NO admissible multiplier structure with <= {MAXO} orbits"))
    print(f"\n{good} entries have some admissible structure, {bad} have none")

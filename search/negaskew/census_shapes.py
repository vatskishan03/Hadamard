"""Plug-in shapes (K = Z_m, q prime = 3 mod 4, gcd(m, q) = 1, H_K odd order) for the census-undecided v (odd v < 2500 without a
matrix of order 4v, search/data/skew_census_undecided.txt, minus the orders of skew/): the shape list of the real-class sweep
sweep_real_census.py.   usage: python3 census_shapes.py       NUMAX (default 10) K-orbits at most; PLUG_SETTLED=repo as in enum_plug.py
Writes runs/census_shapes_n<NUMAX>.json (rows (nu, v, m, q, gens, |H_K|)) and prints the number of shapes per nu.
module: shapes(numax)"""
import json, sys, os
from math import gcd
from sympy import isprime
from enum_plug import census_values, subgroups_odd, orbits
def shapes(numax):
    shapes = []
    for v in census_values():
        for q in range(3, v):
            if v % q or not isprime(q) or q % 4 != 3: continue
            m = v // q
            if m % 2 == 0 or m < 3 or gcd(m, q) != 1: continue
            subs = subgroups_odd(m)
            for H, gens in subs.items():
                nu = len(orbits(m, H))
                if nu < 3 or nu > numax: continue
                shapes.append((nu, v, m, q, list(gens), len(H)))
    shapes.sort()
    return shapes
if __name__ == "__main__":
    NUMAX = int(os.environ.get('NUMAX', '10'))
    sh = shapes(NUMAX)
    os.makedirs('runs', exist_ok=True)
    json.dump(sh, open(f'runs/census_shapes_n{NUMAX}.json', 'w'))
    from collections import Counter
    print(len(census_values()), 'census values not in skew/;', len(sh), 'shapes', sorted(Counter(s[0] for s in sh).items()))

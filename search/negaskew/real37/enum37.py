"""enumerate real plug-in shapes that the sweep of sweep_real_census.py (K = Z_m, q prime, gcd(m, q) = 1) does not cover, for the census values not settled in skew/: K cyclic with q | m, K non-cyclic (diagonal-unit H_K), prime-power q (27, 343, ...).
Output shapes37.json: [nu, v, Kmods, q, gens, Hsize, family] (family A: q prime and q | m, C: q a prime power, B: non-cyclic K).  Environment NUMAX (default 11): the largest number of K-orbits.  usage: python3 enum37.py"""
import json, os, sys, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.join(HERE, '..', '..', '..')
from math import gcd
from sympy import isprime, factorint
from collections import Counter
sys.argv = sys.argv[:1]
NUMAX = int(os.environ.get('NUMAX', '11'))
und = [int(l) for l in open(os.path.join(REPO, 'search', 'data', 'skew_census_undecided.txt')) if l.strip() and not l.startswith('#')]
sk = [int(d) for d in os.listdir(os.path.join(REPO, 'skew')) if d.isdigit()]
settled = {N // 4 for N in sk if N % 8 == 4} | {933, 919}
und2 = [v for v in und if v not in settled]
def odd_subgroups_cyclic(m):
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
    for u in odd: subs.setdefault(closure([u]), [u])
    base = list(subs.items())
    for (H1, g1), (H2, g2) in itertools.combinations(base, 2):
        if H1 <= H2 or H2 <= H1: continue
        H = closure(g1 + g2)
        if len(H) % 2 == 0: continue
        subs.setdefault(H, g1 + g2)
    return subs
def n_orbits_cyclic(m, H):
    seen = [False] * m; n = 0
    for k in range(m):
        if not seen[k]:
            n += 1
            for h in H: seen[k * h % m] = True
    return n
def ppower(q):
    f = factorint(q)
    return len(f) == 1
shapes = []
for v in und2:
    # q: prime or prime power = 3 mod 4 dividing v
    qs = set()
    for p in factorint(v):
        k = 1
        while v % (p ** k) == 0:
            q = p ** k
            if q % 4 == 3 and q >= 3 and q < v: qs.add(q)
            k += 1
    for q in sorted(qs):
        m = v // q
        if m % 2 == 0 or m < 3: continue
        p = list(factorint(q))[0]
        # family A / C: K = Z_m cyclic (gcd(m,q) may be > 1, q may be a prime power)
        fam = 'A' if isprime(q) else 'C'
        if gcd(m, q) == 1 and isprime(q): pass     # K cyclic, q prime, gcd(m, q) = 1: covered by ../sweep_real_census.py
        else:
            for H, gens in odd_subgroups_cyclic(m).items():
                nu = n_orbits_cyclic(m, H)
                if 3 <= nu <= NUMAX: shapes.append((nu, v, [m], q, list(gens), len(H), fam))
        # family B: non-cyclic K of order m: Kmods = prime-power components, grouped: only when some prime appears with exponent >= 2
        f = factorint(m)
        if any(e >= 2 for e in f.values()):
            # abelian groups of order m: for each prime p^e choose a partition of e; cyclic coprime parts merged into one modulus
            from sympy.utilities.iterables import partitions
            per_prime = []
            for pr, e in f.items():
                opts = []
                for part in partitions(e):
                    comps = []
                    for sz, mult in part.items(): comps += [pr ** sz] * mult
                    opts.append(sorted(comps))
                per_prime.append(opts)
            for combo in itertools.product(*per_prime):
                # skip the all-cyclic case (covered above): non-cyclic iff some prime has >= 2 components
                if all(len(c) == 1 for c in combo): continue
                # moduli: merge one component of each prime into a single cyclic modulus; extra components separately
                comps = [list(c) for c in combo]
                main = 1
                for c in comps: main *= c[0]
                Kmods = [main] + [x for c in comps for x in c[1:]]
                Kmods.sort(reverse=True)
                units = [[u for u in range(1, mm) if gcd(u, mm) == 1] or [1] for mm in Kmods]
                tot = 1
                for u in units: tot *= len(u)
                if tot > 200000: continue
                seenH = set()
                for u in itertools.product(*units):
                    # order of u
                    x = list(u); o = 1
                    while any(a != 1 for a in x) and o < 10000:
                        x = [a * b % mm for a, b, mm in zip(x, u, Kmods)]; o += 1
                    if o % 2 == 0 or o == 1: continue
                    # subgroup
                    Hs = {tuple([1] * len(Kmods))}; y = tuple(u)
                    cur = tuple([1] * len(Kmods))
                    G = []
                    while True:
                        cur = tuple(a * b % mm for a, b, mm in zip(cur, u, Kmods)); G.append(cur)
                        if all(a == 1 for a in cur): break
                    key = frozenset(G)
                    if key in seenH: continue
                    seenH.add(key)
                    # count orbits on K
                    elems = list(itertools.product(*[range(mm) for mm in Kmods])); seen = set(); nu = 0
                    for e in elems:
                        if e in seen: continue
                        nu += 1
                        for h in G: seen.add(tuple(a * b % mm for a, b, mm in zip(e, h, Kmods)))
                    if 3 <= nu <= NUMAX: shapes.append((nu, v, Kmods, q, [list(u)], len(G), 'B'))
shapes.sort(key=lambda s: (s[0], s[1]))
json.dump(shapes, open('shapes37.json', 'w'))
print(len(shapes), 'new shapes', Counter(s[0] for s in shapes), Counter(s[6] for s in shapes))
print(len(set(s[1] for s in shapes)), 'distinct v')

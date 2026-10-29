"""K = Z_p^2 with a MATRIX group H_K <= GL(2,p) of odd order (Singer subgroups, diagonal tori, unipotent products): shape file for qp37r in the same format as qshape37 (units = the diagonal units that NORMALISE H_K, group perms = translations by H_K-fixed points + reflections).
Note (search/NOTES.md section 90): diag(u1, u2) acts on the K-orbits only if it normalises H_K, which fails for u1 != u2 and a Singer group; writing all diagonal unit pairs would let uA=1 remove solutions.  Only the normalising units are written (the EMPTY verdicts of 575, 1183, 1083 agree with the unreduced runs).
usage: python3 qshape37m.py p q '[[a,b],[c,d]]' '[[...]]' ... OUT.txt   (generators as 2x2 matrices, acting on column vectors (x, y))"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)      # ../qshape.py is the base shape code
import itertools, math, cmath, sys, json
from math import gcd
from qshape import write_shape as _ws
from qshape37 import unit_perms, group_perms
def kshape_matrix(p, mats):
    Km = [p, p]; elems = list(itertools.product(range(p), range(p))); idx = {e: i for i, e in enumerate(elems)}; Ksz = p * p
    def act(M, e): return ((M[0][0] * e[0] + M[0][1] * e[1]) % p, (M[1][0] * e[0] + M[1][1] * e[1]) % p)
    # group generated (as maps on elements: tuple image of each element)
    ident = tuple(range(Ksz))
    maps = {ident}; fr = [ident]
    gm = [tuple(idx[act(M, e)] for e in elems) for M in mats]
    while fr:
        x = fr.pop()
        for g in gm:
            y = tuple(g[x[i]] for i in range(Ksz))
            if y not in maps: maps.add(y); fr.append(y)
    H = sorted(maps)
    if len(H) % 2 == 0: raise ValueError('even order')
    neg = tuple(idx[((-a) % p, (-b) % p)] for (a, b) in elems)
    if neg in H: raise ValueError('-1 in H_K')
    orb = [-1] * Ksz; reps = []; zero = idx[(0, 0)]
    for i in [zero] + [i for i in range(Ksz) if i != zero]:
        if orb[i] < 0:
            o = len(reps); reps.append(i)
            for h in H: orb[h[i]] = o
    nu = len(reps); ell = [orb.count(o) for o in range(nu)]
    negorb = [orb[neg[reps[o]]] for o in range(nu)]
    N = [[[0] * nu for _ in range(nu)] for _ in range(nu)]
    for j in range(nu):
        r = elems[reps[j]]
        for i in range(Ksz):
            k = elems[i]; kk = ((k[0] + r[0]) % p, (k[1] + r[1]) % p)
            N[j][orb[i]][orb[idx[kk]]] += 1
    C = [[0j] * nu for _ in range(nu)]
    for j in range(nu):
        r = elems[reps[j]]
        for i in range(Ksz):
            k = elems[i]; C[j][orb[i]] += cmath.exp(2j * math.pi * ((r[0] * k[0] + r[1] * k[1]) % p) / p)
    return dict(Kmods=Km, nu=nu, Ksz=Ksz, ell=ell, neg=negorb, N=N, C=C, orb=orb, reps=reps, elems=elems, H=H)
def singer(p, order):
    """generator matrix of the subgroup of order `order` (odd, dividing p^2-1) of the Singer cycle F_{p^2}* acting on Z_p^2 = {a + b x}, x^2 = c (c a non-residue)"""
    c = next(c for c in range(2, p) if pow(c, (p - 1) // 2, p) == p - 1)
    def mul(u, v): return ((u[0] * v[0] + c * u[1] * v[1]) % p, (u[0] * v[1] + u[1] * v[0]) % p)
    def power(u, n):
        r = (1, 0)
        while n:
            if n & 1: r = mul(r, u)
            u = mul(u, u); n >>= 1
        return r
    full = p * p - 1
    for a in range(p):
        for b in range(1, p):
            g = (a, b)
            if all(power(g, full // f) != (1, 0) for f in {x for x in range(2, full + 1) if full % x == 0 and all(x % y for y in range(2, int(x ** .5) + 1))}):
                h = power(g, full // order)
                return [[h[0], (c * h[1]) % p], [h[1], h[0]]]
    raise ValueError
def normalising_unit_perms(sh):
    """permutations of the K-orbits induced by the diagonal units u = (u1, u2) of Z_p^2 that normalise H_K (u h u^-1 in H_K for every h); only these map H_K-orbits to H_K-orbits.
    sh must come from kshape_matrix (sh['H'] = the maps of H_K as tuples of element indices)."""
    Km = sh['Kmods']; elems = sh['elems']; orb = sh['orb']; reps = sh['reps']; nu = sh['nu']
    idx = {e: i for i, e in enumerate(elems)}; Hs = set(sh['H'])
    perms = set()
    lists = [[u for u in range(1, m) if gcd(u, m) == 1] or [1] for m in Km]
    for u in itertools.product(*lists):
        d = tuple(idx[tuple((a * x) % m for a, x, m in zip(e, u, Km))] for e in elems)
        dinv = [0] * len(d)
        for i, j in enumerate(d): dinv[j] = i
        if all(tuple(d[h[dinv[i]]] for i in range(len(d))) in Hs for h in sh['H']):
            perms.add(tuple(orb[d[reps[o]]] for o in range(nu)))
    perms.discard(tuple(range(nu)))
    return sorted(perms)
def write_shape_mat(sh, q, fname):
    _ws(sh, q, fname)
    with open(fname, 'a') as f:
        P = normalising_unit_perms(sh); f.write(f"{len(P)}\n")
        for pp in P: f.write(" ".join(map(str, pp)) + "\n")
        G = group_perms(sh); f.write(f"{len(G)}\n")
        for refl, pp in G: f.write(str(refl) + " " + " ".join(map(str, pp)) + "\n")
if __name__ == '__main__':
    p = int(sys.argv[1]); q = int(sys.argv[2]); mats = [json.loads(a) for a in sys.argv[3:-1]]; fn = sys.argv[-1]
    sh = kshape_matrix(p, mats); write_shape_mat(sh, q, fn); print('nu', sh['nu'], 'ell', sh['ell'], '|H_K|', len(sh['H']))

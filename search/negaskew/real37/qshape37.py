"""qshape37: ../qshape.py + the unit-action permutations of the K-orbits appended to the shape file (for the A-canonicalisation uA=1 of qp37r).
Units: diagonal unit tuples u of K = prod Z_m (gcd(u_i, m_i) = 1) acting on K; they commute with H_K, so they permute the K-orbits; a skew block A -> u.A with all other blocks transformed alike is a symmetry of the K-level system.
File tail: 'nU' then nU lines (perm[a] = orbit of u.rep[a]); identity and duplicates removed."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)      # ../qshape.py is the base shape code
import itertools, sys, os, json
from math import gcd
from qshape import kshape, write_shape as _ws, cyc_gens
def unit_perms(sh):
    Km = sh['Kmods']; elems = sh['elems']; orb = sh['orb']; reps = sh['reps']; nu = sh['nu']
    idx = {e: i for i, e in enumerate(elems)}
    perms = set()
    lists = [[u for u in range(1, m) if gcd(u, m) == 1] or [1] for m in Km]
    cnt = 0
    for u in itertools.product(*lists):
        p = tuple(orb[idx[tuple((a * x) % m for a, x, m in zip(elems[reps[o]], u, Km))]] for o in range(nu))
        perms.add(p)
    ident = tuple(range(nu)); perms.discard(ident)
    return sorted(perms)
def group_perms(sh):
    """permutations of the K-orbits induced by translations by H_K-fixed points k0 (singleton orbits) and by reflection composed with them (k -> -k + k0); identity removed.
    A free block X may be replaced by any of these images (EXACT symmetry of the K-level system); together with negation they form the group used by qp37r tF=1."""
    Km = sh['Kmods']; elems = sh['elems']; orb = sh['orb']; reps = sh['reps']; nu = sh['nu']; ell = sh['ell']
    idx = {e: i for i, e in enumerate(elems)}
    fixpts = [reps[o] for o in range(nu) if ell[o] == 1]
    perms = {}
    for k0i in fixpts:
        k0 = elems[k0i]
        for refl in (0, 1):
            p = []
            for o in range(nu):
                e = elems[reps[o]]
                if refl: e = tuple((-a) % m for a, m in zip(e, Km))
                e2 = tuple((a + b) % m for a, b, m in zip(e, k0, Km))
                p.append(orb[idx[e2]])
            perms[tuple(p)] = refl if tuple(p) not in perms else perms[tuple(p)]
    perms.pop(tuple(range(nu)), None)
    return sorted((refl, p) for p, refl in perms.items())
def write_shape(sh, q, fname, units=True):
    _ws(sh, q, fname)
    if units:
        P = unit_perms(sh)
        with open(fname, 'a') as f:
            f.write(f"{len(P)}\n")
            for p in P: f.write(" ".join(map(str, p)) + "\n")
            G = group_perms(sh)
            f.write(f"{len(G)}\n")
            for refl, p in G: f.write(str(refl) + " " + " ".join(map(str, p)) + "\n")
if __name__ == '__main__':
    Km = json.loads(sys.argv[1]); q = int(sys.argv[2]); hg = json.loads(sys.argv[3]); fn = sys.argv[4]
    sh = kshape(Km, hg); write_shape(sh, q, fn); print('nu', sh['nu'], 'units', len(unit_perms(sh)))

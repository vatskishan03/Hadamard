"""Enumerate non-cyclic abelian groups G of order v (homocyclic Sylow parts (Z/p^b)^k, k<=3) and odd-order
multiplier subgroups H <= Aut(G) (cyclic, products, diagonal cyclic, pair joins); compute orbit count n, estimator, DC test.
usage: python3 ncenum.py NMAX out.jsonl v1 v2 ...      appends one json record per shape with n <= NMAX orbits to out.jsonl (e.g. runs/nc_all.jsonl)
Each record: v, the parts, the automorphisms (part, matrix, order), n, orbit sizes, the Gaussian estimate (est, bits) and whether a row-sum pattern exists (dc).
Over all Table 5 entries m >= 4 / census values with a square factor, every shape with n <= 17 orbits has an estimate <= -4.7 bits or no row-sum pattern
(841 = Z_29^2: 9 orbits, estimate +2.3, no pattern); the parts too large to enumerate are 1029, 1369, 1625, 1681, 1849, 2125, 2197 and 2375 (NOTES.md section 88).
Example: python3 ncenum.py 17 runs/nc_539.jsonl 539    -> Z_7^2 x Z_11 with H of order 21 x 5: n = 15, estimate -12.23, dc true (about 1 s)."""
import sys, json, itertools, math, os, time
import numpy as np
from math import gcd
from gshape import *
from dcnp import plain_dc
from sympy import factorint

def odd_matrices_k2(m):
    ent = np.arange(m)
    a, b, c, d = np.meshgrid(ent, ent, ent, ent, indexing='ij')
    a = a.ravel(); b = b.ravel(); c = c.ravel(); d = d.ravel()
    p = factorint(m); pp = list(p.keys())[0]
    det = (a * d - b * c) % m
    ok = (det % pp) != 0
    a, b, c, d, det = a[ok], b[ok], c[ok], d[ok], det[ok]
    A, B, C, D = a.copy(), b.copy(), c.copy(), d.copy()
    order = np.zeros(len(a), dtype=np.int64)
    done = np.zeros(len(a), bool)
    for o in range(1, 2 * m * m + 2):
        ident = (A == 1) & (B == 0) & (C == 0) & (D == 1)
        newly = ident & ~done
        order[newly] = o; done |= newly
        if done.all(): break
        A2 = (A * a + B * c) % m; B2 = (A * b + B * d) % m
        C2 = (C * a + D * c) % m; D2 = (C * b + D * d) % m
        A, B, C, D = A2, B2, C2, D2
    odd = (order % 2 == 1) & done
    tr = (a + d) % m
    # one rep per (order, trace, det)
    keys = {}
    idx = np.nonzero(odd)[0]
    for i in idx:
        k = (int(order[i]), int(tr[i]), int(det[i]))
        if k not in keys: keys[k] = i
    res = []
    for k, i in keys.items():
        res.append((np.array([[a[i], b[i]], [c[i], d[i]]], dtype=np.int64), k[0]))
    return res

def odd_matrices_k3_small(m):
    # brute force for m = 3 only (3^9 matrices)
    res = {}
    I = np.eye(3, dtype=np.int64)
    for ent in itertools.product(range(m), repeat=9):
        M = np.array(ent, dtype=np.int64).reshape(3, 3)
        det = int(round(np.linalg.det(M))) % m
        if det % m == 0: continue
        X = M.copy(); o = 1
        while not np.array_equal(X, I):
            X = (X @ M) % m; o += 1
        if o % 2 == 1:
            key = (o, int(np.trace(M) % m), det, tuple(np.sort(np.linalg.matrix_power(M, 1).ravel() % m)[:0]))
            if key not in res: res[key] = (M, o)
    return list(res.values())

def part_candidates(p, b, k):
    m = p ** b
    if k == 1:
        out = []
        phi = (p - 1) * p ** (b - 1)
        # generator of Z_m^*: find elements of each odd order
        units = [u for u in range(1, m) if gcd(u, m) == 1]
        seen = set()
        for u in units:
            o = 1; x = u
            while x != 1: x = x * u % m; o += 1
            if o % 2 == 1 and o > 1 and o not in seen:
                seen.add(o); out.append((np.array([[u]]), o))
        return out
    if k == 2:
        if m > 31: return None
        return [(M, o) for M, o in odd_matrices_k2(m) if o > 1]
    if k == 3 and m == 3:
        return [(M, o) for M, o in odd_matrices_k3_small(m) if o > 1]
    return None

def orb_key_part(m, k, M):
    mods = [m] * k
    perm = perm_from_matrix(mods, M)
    orb, reps = orbits_of([perm], m ** k)
    szs = sorted(int((orb == i).sum()) for i in range(len(reps)))
    return tuple(szs)

def dedupe_part(m, k, cands):
    seen = {}; out = []
    for M, o in cands:
        key = (o, orb_key_part(m, k, M))
        if key not in seen:
            seen[key] = 1; out.append((M, o, key))
    return out

def group_types(v):
    f = sorted(factorint(v).items())
    per = []
    for p, a in f:
        opts = []
        for b in range(1, a + 1):
            if a % b == 0:
                k = a // b
                if k <= 3 and (p ** b) ** k <= 2200: opts.append((p, b, k))
        per.append(opts)
    return [c for c in itertools.product(*per)]

def enumerate_v(v, nmax, outf, done):
    t0 = time.time()
    for gt in group_types(v):
        if all(k == 1 for (p, b, k) in gt): continue   # cyclic: covered by sweep.py
        parts = [(p ** b, k) for (p, b, k) in gt]
        gkey = json.dumps(parts)
        # per-part candidate lists
        plist = []; ok = True
        for (p, b, k) in gt:
            c = part_candidates(p, b, k)
            if c is None: ok = False; break
            plist.append(dedupe_part(p ** b, k, c))
        if not ok:
            rec = dict(v=v, parts=parts, status='skipped_part_too_large'); open(outf, 'a').write(json.dumps(rec) + '\n'); continue
        mods = []
        for (m, k) in parts: mods += [m] * k
        V = v
        # embedded single generators: list of (descr, perm)
        gens = []
        for j, cl in enumerate(plist):
            for M, o, key in cl:
                mats = [None] * len(parts); mats[j] = M
                _, perm = block_diag_perm(parts, mats)
                gens.append(((j, M.tolist(), o), perm))
        cands = []  # (description, [perms])
        for g in gens: cands.append(([g[0]], [g[1]]))
        # diagonal cyclic: products of gens from different parts
        for (d1, p1), (d2, p2) in itertools.combinations(gens, 2):
            if d1[0] == d2[0]:
                # same part: pair join only if commute
                if np.array_equal(p1[p2], p2[p1]): cands.append(([d1, d2], [p1, p2]))
                continue
            cands.append(([d1, d2], [p1, p2]))          # product group
            if gcd(d1[2], d2[2]) > 1: cands.append(([d1, d2, 'diag'], [p1[p2]]))   # diagonal cyclic
        seen = set()
        cnt = 0
        for desc, perms in cands:
            orb, reps = orbits_of(perms, V); n = len(reps)
            if n > nmax: continue
            szs = tuple(sorted(int((orb == i).sum()) for i in range(n)))
            # require odd group order: check generated group size when small
            key = (n, szs)
            d = shape_g(mods, perms)
            e = est_d(d)
            k2 = (n, szs, round(e, 3))
            if k2 in seen: continue
            seen.add(k2)
            dc = plain_dc(V, d['sz'], d['neg'], n)
            rec = dict(v=v, parts=parts, mods=mods, desc=[(x if x == 'diag' else (x[0], x[1], x[2])) for x in desc], n=n, sz=sorted(set(szs)), est=round(e, 2), dc=bool(dc))
            # store the generator perms' matrices only (descr); perms are rebuilt later
            open(outf, 'a').write(json.dumps(rec) + '\n'); cnt += 1
        print(v, parts, 'candidates', len(cands), 'shapes<=nmax', cnt, 'time %.0f' % (time.time() - t0), flush=True)

if __name__ == "__main__":
    nmax = int(sys.argv[1]); outf = sys.argv[2]
    os.makedirs(os.path.dirname(outf) or '.', exist_ok=True)
    for v in sys.argv[3:]:
        enumerate_v(int(v), nmax, outf, set())

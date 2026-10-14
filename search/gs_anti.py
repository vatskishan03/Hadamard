"""Goethals-Seidel skew-Hadamard matrices of order 4n with an involutive anti-automorphism P (P S P = S^T, P a
symmetric signed permutation), i.e. symmetric orthogonal designs x P + y (S - I) P of type (1, 4n - 1).

S = [[A, BR, CR, DR], [-BR, A, D^T R, -C^T R], [-CR, -D^T R, A, B^T R], [-DR, C^T R, -B^T R, A]] (verify.goethals_seidel)
with A of skew type. P is block monomial: P[i][pi(i)] = eps_i R_{c_i}, R_c the reflection g -> c - g (the diagonal
condition Q A Q^T = A^T forces reflections), pi an involution of the four blocks, eps_{pi(i)} = eps_i,
c_{pi(i)} = c_i. Every entry of P S P = S^T is an equation +-x = +-y between sequence entries; a signed union-find
gives the induced structure of (A, B, C, D) or a contradiction. Block translations diag(T^{a_i}) with sum a_i = 0
keep the Goethals-Seidel form and move c_i to c_i + a_i + a_{pi(i)}, so the centres are normalised.

Usage: python3 search/gs_anti.py n            (lists the consistent structures and their free parameters)
"""
import itertools
import sys


def structures(n):
    # variables: A(g) for g in Z_n (index 0..n-1), B: n.., C: 2n.., D: 3n..
    def var(seq, g):
        return seq * n + g % n
    # entry of S at (block i, row r; block j, col s) as (sign, var) or constant
    def S(i, r, j, s):
        if i == j:
            return (1, var(0, s - r))                    # A circulant: A[r][s] = a(s - r)
        # back-circulant X R: (X R)[r][s] = X[r][n-1-s]... use R = reflection g -> -g: (X R)[r][s] = x(-s - r)
        g = (-s - r) % n
        table = {(0, 1): (1, 1, False), (0, 2): (1, 2, False), (0, 3): (1, 3, False),
                 (1, 0): (-1, 1, False), (1, 2): (1, 3, True), (1, 3): (-1, 2, True),
                 (2, 0): (-1, 2, False), (2, 1): (-1, 3, True), (2, 3): (1, 1, True),
                 (3, 0): (-1, 3, False), (3, 1): (1, 2, True), (3, 2): (-1, 1, True)}
        sg, seq, tr = table[(i, j)]
        if tr:                                           # X^T R: (X^T)[r][s'] = x(r - s')
            g = (-g) % n
        return (sg, var(seq, g))

    out = []
    invols = [p for p in itertools.permutations(range(4)) if all(p[p[i]] == i for i in range(4))]
    for pi in invols:
        orbits = sorted({tuple(sorted((i, pi[i]))) for i in range(4)})
        for eps_o in itertools.product((1, -1), repeat=len(orbits)):
            if eps_o[0] == -1:
                continue                                  # global sign of P is irrelevant
            eps = [0] * 4
            for o, e in zip(orbits, eps_o):
                for i in o:
                    eps[i] = e
            # centres: normalise as described; enumerate the remaining freedom (at most one centre free)
            for cs in itertools.product(range(n), repeat=len(orbits)):
                c = [0] * 4
                for o, x in zip(orbits, cs):
                    for i in o:
                        c[i] = x
                if not normalised(pi, c, n):
                    continue
                res = check(n, pi, eps, c, S)
                if res is not None:
                    out.append((pi, tuple(eps), tuple(c), res))
    return out


def normalised(pi, c, n):
    # block translations with sum 0: c_i -> c_i + a_i + a_pi(i). Fixed blocks: c_i can be made 0 (2 invertible),
    # transposition (i j): c_i + a_i + a_j. Keep only centre vectors that are lexicographically minimal in their class
    # under these moves (brute force over a: n^3 choices is too many; use the invariant instead).
    fixed = [i for i in range(4) if pi[i] == i]
    pairs = sorted({tuple(sorted((i, pi[i]))) for i in range(4) if pi[i] != i})
    # invariant: the moves shift fixed centres by 2 a_i and pair centres by a_i + a_j, subject to sum a = 0.
    # Sum over all blocks of the shift contributions: sum_fixed 2 a_i + sum_pairs 2 (a_i + a_j) = 2 sum a = 0, so
    # I = sum_fixed c_i + 2 sum_pairs c_i (each pair counted once, doubled) is invariant (mod n), and every
    # configuration with the same I is reachable. Normal form: all centres 0 except the last orbit, set to make I.
    orbits = sorted({tuple(sorted((i, pi[i]))) for i in range(4)})
    for o in orbits[:-1]:
        if c[o[0]] != 0:
            return False
    return True


def check(n, pi, eps, c, S):
    N = 4 * n
    parent = list(range(4 * n))
    par = [0] * (4 * n)                                   # parity to parent: x = (-1)^par * parent

    def find(x):
        if parent[x] == x:
            return x, 0
        r, p = find(parent[x])
        parent[x] = r
        par[x] ^= p
        return r, par[x]

    def union(x, y, s):                                  # x = s * y
        rx, px = find(x)
        ry, py = find(y)
        want = 0 if s == 1 else 1
        if rx == ry:
            return (px ^ py) == want
        parent[rx] = ry
        par[rx] = px ^ py ^ want
        return True

    # skew type of A: a(0) = 1 (constant), a(-g) = -a(g)
    for g in range(1, n):
        if not union(g, (-g) % n, -1):
            return None
    sys.setrecursionlimit(10000)
    for i in range(4):
        for r in range(n):
            pr = (c[i] - r) % n                            # P maps row (i, r) to (pi(i), c_i - r)
            for j in range(4):
                for s in range(n):
                    ps = (c[j] - s) % n
                    lhs_sign = eps[i] * eps[j]
                    e1 = S(pi[i], pr, pi[j], ps)
                    e2 = S(j, s, i, r)                    # S^T[(i,r),(j,s)] = S[(j,s),(i,r)]
                    s1, v1 = e1
                    s2, v2 = e2
                    s1 *= lhs_sign
                    if v1 == 0 and v2 == 0:              # a(0) = 1 both sides
                        if s1 != s2:
                            return None
                        continue
                    if v1 == 0 or v2 == 0:
                        # a(0) equated to another variable: fine, it fixes that variable
                        x, y = (v2, v1) if v1 == 0 else (v1, v2)
                        sgn = s1 * s2
                        if not union(x, y, sgn):
                            return None
                        continue
                    if not union(v1, v2, s1 * s2):
                        return None
    roots = {}
    for x in range(4 * n):
        r, p = find(x)
        roots.setdefault(r, []).append((x, p))
    root0 = find(0)[0]
    free = [r for r in roots if r != root0]
    return len(free), roots


def main():
    n = int(sys.argv[1])
    res = structures(n)
    print(f'n = {n}: {len(res)} consistent (pi, eps, c) configurations')
    seen = {}
    for pi, eps, c, (nfree, roots) in res:
        key = (pi, eps, nfree)
        seen.setdefault(key, []).append(c)
    for (pi, eps, nfree), cs in sorted(seen.items(), key=lambda x: x[0][2]):
        print(f'  pi={pi} eps={eps} free classes {nfree}  centres {cs[:4]}{"..." if len(cs) > 4 else ""} ({len(cs)})')


if __name__ == '__main__':
    main()


def search(n, maxfree=26, verbose=True):
    """Every consistent configuration: enumerate the free classes, test sum PAF = 4n delta."""
    import numpy as np
    res = structures(n)
    found = []
    for pi, eps, c, (nfree, roots) in res:
        if nfree > maxfree:
            continue
        root0 = None
        # identify the class of a(0) (fixed +1)
        for r, members in roots.items():
            if any(x == 0 for x, p in members):
                root0 = r
        free = [r for r in roots if r != root0]
        # value of each variable = (-1)^p * value(root)
        idx = np.zeros(4 * n, dtype=np.int64)
        sgn = np.zeros(4 * n, dtype=np.int64)
        for k, r in enumerate(free):
            for x, p in roots[r]:
                idx[x] = k + 1
                sgn[x] = -1 if p else 1
        for x, p in roots[root0]:
            idx[x] = 0
            sgn[x] = -1 if p else 1
        # sign of a(0)'s class: a(0) = +1 means root value = sgn[0]
        base = sgn[0]
        nf = len(free)
        cnt = 0
        for bits in range(1 << nf):
            vals = np.empty(nf + 1, dtype=np.int64)
            vals[0] = base
            for k in range(nf):
                vals[k + 1] = -1 if (bits >> k) & 1 else 1
            x = sgn * vals[idx]
            seqs = x.reshape(4, n)
            tot = np.zeros(n, dtype=np.int64)
            for q in seqs:
                tot += np.array([q @ np.roll(q, -s) for s in range(n)])
            if tot[0] == 4 * n and not tot[1:].any():
                cnt += 1
                found.append((pi, eps, c, seqs.copy()))
        if verbose and cnt:
            print(f'  pi={pi} eps={eps} c={c}: {cnt} GS quadruples with this anti-automorphism ({nf} free)')
    return found


if __name__ == '__main__' and len(sys.argv) > 2 and sys.argv[2] == 'search':
    f = search(int(sys.argv[1]))
    print('total', len(f))


def class_layout(n, roots):
    """For each free class: the set of sequences (0 A, 1 B, 2 C, 3 D) its members lie in."""
    root0 = next(r for r, mem in roots.items() if any(x == 0 for x, p in mem))
    lay = []
    for r, mem in roots.items():
        if r == root0:
            continue
        lay.append(sorted({x // n for x, p in mem}))
    return lay


def search_fast(n, maxfree=24, chunk=1 << 16):
    """Vectorised version of search() for configurations with at most maxfree free classes."""
    import numpy as np
    res = structures(n)
    found, big = [], []
    for pi, eps, c, (nfree, roots) in res:
        if nfree > maxfree:
            big.append((pi, eps, c, nfree, roots))
            continue
        root0 = next(r for r, mem in roots.items() if any(x == 0 for x, p in mem))
        free = [r for r in roots if r != root0]
        idx = np.zeros(4 * n, dtype=np.int64)
        sgn = np.zeros(4 * n, dtype=np.int64)
        for k, r in enumerate(free):
            for x, p in roots[r]:
                idx[x], sgn[x] = k + 1, (-1 if p else 1)
        for x, p in roots[root0]:
            idx[x], sgn[x] = 0, (-1 if p else 1)
        base = sgn[0]
        nf = len(free)
        for start in range(0, 1 << nf, chunk):
            b = np.arange(start, min(1 << nf, start + chunk), dtype=np.int64)
            vals = np.empty((len(b), nf + 1), dtype=np.int64)
            vals[:, 0] = base
            for k in range(nf):
                vals[:, k + 1] = 1 - 2 * ((b >> k) & 1)
            X = (sgn[None, :] * vals[:, idx]).reshape(len(b), 4, n)
            ok = np.ones(len(b), dtype=bool)
            for s in range(1, n // 2 + 1):
                tot = (X * np.roll(X, -s, axis=2)).sum(axis=(1, 2))
                ok &= tot == 0
                if not ok.any():
                    break
            for k in np.flatnonzero(ok):
                found.append((pi, eps, c, X[k].copy()))
    return found, big

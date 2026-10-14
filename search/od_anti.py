"""Involutive anti-automorphisms of the ingredients of Seberry's Theorem 7 (skew-amicable-plug certificates).

A symmetric signed permutation P1 with P1 (aI + bB + cC) P1 = aI - bB - cC turns the skew-type OD of the theorem into
the symmetric OD a P1 + b P1 B + c P1 C, and then P1 B (x) J + (P1 + P1 C) (x) D is a symmetric Hadamard matrix of
order 4t (h - 1) whenever 4t = m h (the weighted plug-in of NOTES §32; D the core of any symmetric Hadamard matrix
of order h). For skew/7956 that would be a symmetric matrix of order 7956.

  python3 search/od_anti.py pair CERT        signed permutations s of the amicable pair (S, M) of a certificate with
                                             s S s^T = S^T and s M s^T = +-M (backtracking)
  python3 search/od_anti.py cert CERT        block-monomial reflections P1 (blocks eps_i R_{c_i}, any centres) with
                                             P1 G P1 = G^T for the certificate's OD in code form G
  python3 search/od_anti.py structures t     every block-monomial P1 for a Goethals-Seidel OD(4t; 1, 3, 4t - 4) in the
                                             normal form of od13.c (b entries at 0 of blocks 2-4, any signs): the
                                             induced ties between the sequences and the parity test of the bc condition
"""
import itertools, json, sys
import numpy as np
from pathlib import Path as _P
sys.path.insert(0, str(_P(__file__).resolve().parent.parent / 'verify'))
import verify as V


def pair_search(cert):
    n = cert['n']; m = 4 * n
    P, N, _ = V._skew_propus_od(n, cert['blocks'], cert['rowsums'])
    I = np.eye(m, dtype=np.int64)
    S = (I + P @ N).tolist(); M = (P + N).tolist()
    ST = [list(r) for r in zip(*S)]
    neg = lambda X: [[-x for x in r] for r in X]

    def search(pairs):
        sols, pi, s, used = [], [-1] * m, [0] * m, [False] * m

        def rec(k):
            if k == m:
                sols.append((pi[:], s[:])); return
            for t in range(m):
                if used[t] or any(X[k][k] != Y[t][t] for X, Y in pairs):
                    continue
                for si in ((1,) if k == 0 else (1, -1)):
                    if all(si * s[j] * X[k][j] == Y[t][pi[j]] and si * s[j] * X[j][k] == Y[pi[j]][t]
                           for j in range(k) for X, Y in pairs):
                        pi[k], s[k], used[t] = t, si, True
                        rec(k + 1)
                        pi[k], s[k], used[t] = -1, 0, False
        rec(0)
        return sols
    for name, pairs in [('s S s^T = S^T, s M s^T = M', [(S, ST), (M, M)]),
                        ('s S s^T = S^T, s M s^T = -M', [(S, ST), (M, neg(M))]),
                        ('s S s^T = S, s M s^T = M', [(S, S), (M, M)])]:
        print(f'{name}: {len(search(pairs))} signed permutations (up to a global sign)')


def code_matrix(blocks):
    return V.goethals_seidel(*[V.circulant(np.array(b, dtype=np.int64)) for b in blocks])


def invols():
    return [p for p in itertools.permutations(range(4)) if all(p[p[i]] == i for i in range(4))]


def sign_patterns(pi):
    orbits = sorted({tuple(sorted((i, pi[i]))) for i in range(4)})
    for eps_o in itertools.product((1, -1), repeat=len(orbits)):
        if eps_o[0] == -1:
            continue                                         # a global sign of P1 changes nothing
        eps = [0] * 4
        for o, e in zip(orbits, eps_o):
            for i in o:
                eps[i] = e
        yield orbits, tuple(eps)


def cert_search(cert):
    """P1[i][pi(i)] = eps_i Refl_{c_i}, Refl_c: r -> c - r; the condition splits into one table per block pair."""
    blocks = cert['od']['blocks']
    n = len(blocks[0])
    G = code_matrix(blocks)
    blk = lambda a, b: G[a * n:(a + 1) * n, b * n:(b + 1) * n]
    r = np.arange(n)
    ci, cj = np.arange(n)[:, None, None, None], np.arange(n)[None, :, None, None]
    rr, ss = r[None, None, :, None], r[None, None, None, :]
    tab = {}
    for a, b in itertools.product(range(4), repeat=2):
        src = blk(a, b)[(ci - rr) % n, (cj - ss) % n]
        for j, i in itertools.product(range(4), repeat=2):
            tgt = blk(j, i).T[None, None, :, :]
            for sg in (1, -1):
                tab[(a, b, j, i, sg)] = (sg * src == tgt).all(axis=(2, 3))
    found = 0
    for pi in invols():
        for orbits, eps in sign_patterns(pi):
            k = len(orbits)
            oof = {i: next(q for q, o in enumerate(orbits) if i in o) for i in range(4)}
            ok = np.ones((n,) * k, dtype=bool)                  # over the centres of the block orbits
            for i, j in itertools.product(range(4), repeat=2):
                tb = tab[(pi[i], pi[j], j, i, eps[i] * eps[j])]
                a, b = oof[i], oof[j]
                if a == b:
                    m2 = np.diagonal(tb).copy()
                    shape = [1] * k; shape[a] = n
                    ok &= m2.reshape(shape)
                else:
                    shape = [1] * k; shape[a] = n; shape[b] = n
                    ok &= (tb if a < b else tb.T).reshape(shape)
            for cs in zip(*np.nonzero(ok)):
                c = [int(cs[oof[i]]) for i in range(4)]
                found += 1
                print('  P1:', pi, eps, c)
    print(f'block-monomial reflections with P1 G P1 = G^T: {found}')


def structure(t, pi, eps, bsigns, G):
    """Signed union-find of P1 G P1 = G^T on the variables W_k(g) (k = 0..3, g in Z_t) of the od13 normal form, P1
    with blocks eps_i R (R the back identity, forced by the b part). None: contradiction; a string: type clash."""
    R = V.back_identity(t)
    P = np.zeros((4 * t, 4 * t), dtype=np.int64)
    for i in range(4):
        P[i * t:(i + 1) * t, pi[i] * t:(pi[i] + 1) * t] = eps[i] * R
    L = P @ G @ P
    parent, par = list(range(4 * t)), [0] * (4 * t)

    def find(x):
        p = 0
        while parent[x] != x:
            p ^= par[x]; x = parent[x]
        return x, p

    def union(x, y, s):
        (rx, px), (ry, py) = find(x), find(y)
        want = 0 if s == 1 else 1
        if rx == ry:
            return (px ^ py) == want
        parent[rx], par[rx] = ry, px ^ py ^ want
        return True
    for g in range(1, t):                                    # W_0 is skew
        if not union(g, (-g) % t, -1):
            return None
    for u, w in zip(L.ravel(), G.T.ravel()):
        if not union(abs(int(u)) - 1, abs(int(w)) - 1, int(np.sign(u) * np.sign(w))):
            return None
    roots = {}
    for x in range(4 * t):
        rt, p = find(x)
        roots.setdefault(rt, []).append((x, p))
    typ = lambda x: 'a' if x == 0 else ('b' if x % t == 0 else 'c')
    for mem in roots.values():
        ts = {typ(x) for x, p in mem}
        if len(ts) > 1:
            return 'type clash'
        if ts == {'a'} and any(p for x, p in mem):
            return 'a sign'
        if ts == {'b'} and len({p ^ (bsigns[x // t - 1] < 0) for x, p in mem}) > 1:
            return 'b signs'
    return [mem for mem in roots.values() if typ(mem[0][0]) == 'c']


def structures(t):
    G = V.goethals_seidel(*[V.circulant(np.array([k * t + g + 1 for g in range(t)], dtype=np.int64))
                            for k in range(4)])
    counts, feasible = {}, 0
    for bs in itertools.product((1, -1), repeat=3):
        for pi in invols():
            for _, eps in sign_patterns(pi):
                cl = structure(t, pi, eps, bs, G)
                if cl is None or isinstance(cl, str):
                    continue
                desc = set()
                for mem in cl:
                    ks = sorted({x // t for x, p in mem})
                    if ks == [0]:
                        continue
                    (x1, p1), (x2, p2) = mem
                    if len(ks) == 1:
                        desc.add(f'W{ks[0] + 1} {"symmetric" if p1 == p2 else "skew"}')
                    else:
                        g1, g2 = x1 % t, x2 % t
                        rel = 'g' if g1 == g2 else '-g' if (g1 + g2) % t == 0 else f'{(g1 + g2) % t}-g'
                        desc.add(f'W{ks[1] + 1}(g) = {"+" if p1 == p2 else "-"}W{ks[0] + 1}({rel})')
                key = tuple(sorted(desc))
                counts[key] = counts.get(key, 0) + 1
                # bc condition: on every pair {s, -s} the symmetric parts of W2, W3, W4 sum to 0
                sym = [d for d in key if d.endswith('symmetric')]
                ties = [d for d in key if '=' in d]
                if not (len(sym) == 1 and len(ties) == 1):
                    feasible += 1
    print(f't = {t}: {sum(counts.values())} consistent (pi, eps, b signs) configurations')
    for key, c in sorted(counts.items()):
        print(f'  {c:3d} x  ' + '; '.join(key))
    print(f'configurations not of the form "one symmetric block, the other two tied": {feasible}')
    print('In that form the symmetric parts S of the tied pair satisfy S_l = +-S_k, so S_k + S_l is 0 or +-2 and can'
          ' never cancel the +-1 of the symmetric block: the bc condition fails on every pair.')


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'pair':
        pair_search(json.load(open(sys.argv[2])))
    elif cmd == 'cert':
        cert_search(json.load(open(sys.argv[2])))
    elif cmd == 'structures':
        structures(int(sys.argv[2]))

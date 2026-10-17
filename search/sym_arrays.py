#!/usr/bin/env python3
"""Symmetric 4 x 4 arrays of circulant blocks with the reversal R (search/NOTES.md section 49).

Four circulant variables X_1..X_4, each either symmetric (blocks +-X or +-XR) or arbitrary (blocks +-XR or
+-RX = +-X^T R; XR and RX are symmetric for every circulant X, but they do not commute). An array is symmetric iff its
symbol pattern is, and H^2 has the diagonal blocks sum X X^T iff every row holds each variable once, so the pattern is
a symmetric Latin square. Writing every block as M R^e (M = X or X^T), a product is
(M1 R^e1)(M2 R^e2) = M1 M2' R^(e1 + e2) with M2' = M2^T when e1 = 1, and for generic circulants the off-diagonal
blocks of H^2 vanish identically iff the signed monomials cancel. The search assigns blocks row by row and checks
each off-diagonal block of H^2 as soon as both rows are complete.

    python3 search/sym_arrays.py           every type assignment, four distinct variables; then X_3 = X_2 (control)
    python3 search/sym_arrays.py --klein   blocks twisted by any involutive multiplier in {1, -1, u, -u}

Result: no array with four distinct variables, for any choice of which are symmetric; with X_3 = X_2 (Balonin's
propus/Balonin array among them) 768 arrays when X_1 is symmetric and X_2, X_4 arbitrary. With the Klein twists,
still none with four distinct variables (35 type combinations); the twisted repeated-block arrays impose the same
autocorrelation equation as propus.
"""
import itertools
from collections import Counter

def latin_squares():
    out = []
    for rows in itertools.permutations(itertools.permutations(range(4)), 4):
        if all(len({r[j] for r in rows}) == 4 for j in range(4)) and \
           all(rows[i][j] == rows[j][i] for i in range(4) for j in range(4)):
            out.append(rows)
    return out

def count(L, name, symset):
    """Number of sign/orientation choices making the array orthogonal; name maps pattern symbol -> variable."""
    def blocks(v):
        return [(v, 0, 0), (v, 0, 1)] if v in symset else [(v, 0, 1), (v, 1, 1)]
    def prod(b1, b2):
        (v1, t1, e1), (v2, t2, e2) = b1, b2
        t1p = 0 if v1 in symset else t1
        t2p = 0 if v2 in symset else t2 ^ e1
        return (tuple(sorted([(v1, t1p), (v2, t2p)])), (e1 + e2) % 2)
    order = [(i, j) for i in range(4) for j in range(i, 4)]          # row by row
    ent = {}
    def check(i, k):
        c = Counter()
        for j in range(4):
            b1, s1 = ent[(min(i, j), max(i, j))]
            b2, s2 = ent[(min(j, k), max(j, k))]
            c[prod(b1, b2)] += s1 * s2
        return not any(c.values())
    total = 0
    def rec(p):
        nonlocal total
        if p == len(order):
            total += 1
            return
        i, j = order[p]
        for b in blocks(name[L[i][j]]):
            for s in ((1,) if p == 0 else (1, -1)):                 # global sign
                ent[(i, j)] = (b, s)
                # after the last entry of row r (j == 3) every block (i', r) with i' < r is complete
                ok = True
                if j == 3:
                    r = i
                    for i2 in range(r):
                        if not check(i2, r):
                            ok = False
                            break
                if ok:
                    rec(p + 1)
        ent.pop((i, j), None)
    rec(0)
    return total

# Twisted blocks: the Klein group {e, a, b, ab} of involutive multipliers of Z_v (a = -1, b = u, ab = -u, for v with
# two odd prime factors). Variable k has a type t_k and is invariant under the multiplier a*t_k; its blocks are
# X P_t and P_t X = X^T P_t (both symmetric), and (M1 P_t1)(M2 P_t2) = M1 M2^(t1) P_(t1 t2).
# encode K as bit pairs: e=0, a=1, b=2, ab=3 ; product = xor
E, A_, B_, AB = 0, 1, 2, 3
def _unused_latin():
    out = []
    for rows in itertools.permutations(itertools.permutations(range(4)), 4):
        if all(len({r[j] for r in rows}) == 4 for j in range(4)) and all(rows[i][j] == rows[j][i] for i in range(4) for j in range(4)):
            out.append(rows)
    return out
def count_klein(L, name, types):
    # canonical form of X_k^(w) given transposition flag tr (0: X, 1: X^T): returns flag
    def tr_after(k, flag, w):
        t = types[k]; m = A_ ^ t
        if t == E: return 0                      # symmetric: X = X^T
        # X^T = X^(a); so (X^tr)^(w) = X^(a^tr * w)
        ww = w ^ (A_ if flag else 0)
        return 0 if ww in (E, m) else 1
    def norm(k, flag):
        return 0 if types[k] == E else flag
    def blocks(k):
        return [(k, 0), (k, 1)] if types[k] != E else [(k, 0)]
    def prod(b1, b2):
        (k1, f1), (k2, f2) = b1, b2
        t1, t2 = types[k1], types[k2]
        m1 = (k1, norm(k1, f1))
        m2 = (k2, tr_after(k2, f2, t1))
        return (tuple(sorted([m1, m2])), t1 ^ t2)
    order = [(i, j) for i in range(4) for j in range(i, 4)]
    ent = {}
    def check(i, k):
        c = Counter()
        for j in range(4):
            b1, s1 = ent[(min(i, j), max(i, j))]; b2, s2 = ent[(min(j, k), max(j, k))]
            c[prod(b1, b2)] += s1 * s2
        return not any(c.values())
    total = 0; example = None
    def rec(p):
        nonlocal total, example
        if p == len(order):
            total += 1
            if example is None: example = dict(ent)
            return
        i, j = order[p]
        for b in blocks(name[L[i][j]]):
            for s in ((1,) if p == 0 else (1, -1)):
                ent[(i, j)] = (b, s)
                ok = True
                if j == 3:
                    for i2 in range(i):
                        if not check(i2, i): ok = False; break
                if ok: rec(p + 1)
        ent.pop((i, j), None)
    rec(0)
    return total, example

if __name__ == '__main__':
    lat = latin_squares()
    print('symmetric Latin squares of order 4:', len(lat))
    for nsym in range(5):
        tot = 0
        for symset in itertools.combinations(range(4), nsym):
            tot += sum(count(L, [0, 1, 2, 3], set(symset)) for L in lat)
        print(f'four distinct variables, {nsym} of them symmetric: {tot} arrays', flush=True)
    import sys
    if '--klein' in sys.argv:
        for types in itertools.product(range(4), repeat=4):
            if list(types) == sorted(types):
                tot = sum(count_klein(L, [0, 1, 2, 3], types)[0] for L in lat)
                print('Klein types', types, ':', tot, 'arrays with four distinct variables', flush=True)
        print('control X_3 = X_2, types (e, a, a, a):', sum(count_klein(L, [0, 1, 1, 3], (0, 1, 1, 1))[0] for L in lat))
        sys.exit()
    for symset in [{0}, {0, 3}, set(), {0, 1, 3}]:
        tot = sum(count(L, [0, 1, 1, 3], symset) for L in lat)
        print(f'X_3 = X_2, symmetric variables {sorted(symset)}: {tot} arrays', flush=True)

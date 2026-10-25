"""K_{9,9} placements in the Paley two-graph of Z_37 (or any conference core given as a sign matrix) whose sign matrix is
switching-equivalent to H' (x) J_3, H' = [[1,1,1],[1,1,-1],[1,-1,1]] (the only BH(9,6) sign class within Lindsey's
bound): A = A1 u A2 u A3, B = B1 u B2 u B3 (triples), S[a, b] = u_a v_b H'[i(a), j(b)].
Search: A1 = {a0, a1, a2}; B lies in {b : S[a0,b] S[a1,b] = e1, S[a0,b] S[a2,b] = e2}."""
import itertools, sys, pickle
import numpy as np
p = 37
qr = {x * x % p for x in range(1, p)}
S = np.array([[0 if i == j else (1 if (j - i) % p in qr else -1) for j in range(p)] for i in range(p)])
Hp = np.array([[1, 1, 1], [1, 1, -1], [1, -1, 1]])
found = set()
def check(A, B):
    """A, B disjoint 9-sets; test S[A,B] ~ H' (x) J3"""
    M = S[np.ix_(A, B)]
    # normalise: switch rows so column 0 is +, columns so row 0 is +
    M = M * M[:, :1]; M = M * M[:1, :]
    rows = {}
    for i, r in enumerate(M): rows.setdefault(tuple(r), []).append(i)
    if len(rows) != 3 or any(len(v) != 3 for v in rows.values()): return False
    cols = {}
    for j, c in enumerate(M.T): cols.setdefault(tuple(c), []).append(j)
    if len(cols) != 3 or any(len(v) != 3 for v in cols.values()): return False
    R = np.array([list(r) for r in rows])[:, [v[0] for v in cols.values()]]
    G = R @ R.T
    return all(abs(G[i, j]) == 1 for i in range(3) for j in range(i + 1, 3))
others = list(range(1, p))
for a1, a2 in itertools.combinations(others, 2):
    a0 = 0
    for e1 in (1, -1):
        for e2 in (1, -1):
            C = [b for b in range(p) if b not in (a0, a1, a2) and S[a0, b] * S[a1, b] == e1 and S[a0, b] * S[a2, b] == e2]
            if len(C) < 9: continue
            for B in itertools.combinations(C, 9):
                # group the other points by their B-vector up to sign
                groups = {}
                for x in range(p):
                    if x in B or x in (a0, a1, a2): continue
                    v = S[x, list(B)]; v = tuple(v * v[0])
                    groups.setdefault(v, []).append(x)
                vA1 = tuple(S[a0, list(B)] * S[a0, B[0]])
                cand = [g for v, g in groups.items() if len(g) >= 3 and v != vA1]
                for g2, g3 in itertools.combinations(cand, 2):
                    for A2 in itertools.combinations(g2, 3):
                        for A3 in itertools.combinations(g3, 3):
                            A = [a0, a1, a2] + list(A2) + list(A3)
                            if check(A, list(B)):
                                found.add((tuple(sorted(A)), tuple(sorted(B))))
print("placements with 0 in A:", len(found))
allp = set()
for A, B in found:
    for t in range(p):
        allp.add((tuple(sorted((x + t) % p for x in A)), tuple(sorted((y + t) % p for y in B))))
print("all placements (A, B ordered):", len(allp))
pickle.dump(sorted(allp), open("k99_37.pkl", "wb"))

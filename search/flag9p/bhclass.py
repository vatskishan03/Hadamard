"""Sign classes (row/column switching and permutations) of normalised BH(n,6) matrices from bh output lines."""
import sys, itertools
import numpy as np
from pynauty import Graph, certificate
n = int(sys.argv[1]); path = sys.argv[2]
def signs(line):
    k = [0]*n + [int(c) for c in line.strip()]
    M = np.zeros((n, n), int); M[0, :] = 1
    for a in range(1, n):
        M[a, 0] = 1
        for b in range(1, n): M[a, b] = 1 if int(line[(a-1)*(n-1) + b-1]) % 2 == 0 else -1
    return M
def cert(M):
    best = None
    for r0 in range(n):
        for c0 in range(n):
            N = M * M[:, [c0]] * M[[r0], :] * M[r0, c0]   # normalise row r0 and column c0 to +
            rows = [r for r in range(n) if r != r0]; cols = [c for c in range(n) if c != c0]
            adj = {i: [] for i in range(2 * (n - 1))}
            for i, r in enumerate(rows):
                for j, c in enumerate(cols):
                    if N[r, c] > 0: adj[i].append(n - 1 + j)
            g = Graph(2 * (n - 1), directed=False, adjacency_dict=adj, vertex_coloring=[set(range(n - 1)), set(range(n - 1, 2 * (n - 1)))])
            ce = certificate(g)
            if best is None or ce < best: best = ce
    return best
seen = set(); classes = {}
for line in open(path):
    line = line.strip()
    if len(line) != (n - 1) ** 2: continue
    M = signs(line); key = M.tobytes()
    if key in seen: continue
    seen.add(key)
    c = cert(M); classes.setdefault(c, M)
print(f"n={n}: distinct normalised sign matrices {len(seen)}, sign classes {len(classes)}")
import pickle; pickle.dump(list(classes.values()), open(f'bh{n}_signclasses.pkl', 'wb'))
for M in list(classes.values())[:10]:
    print((M > 0).astype(int).tolist())

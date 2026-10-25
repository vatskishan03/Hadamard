"""Sign classes of BH(9, 6): 9x9 +-1 matrices M with all row and all column inner products in {+-3, +-9} (vanishing sums
of nine sixth roots are three triangles or one triangle and three antipodal pairs) that admit cube-root phases.
Normalised: first row and first column +1, so rows have 0, 3 or 6 minus signs among the last 8 coordinates and pairwise
Hamming distances 0, 3 or 6. Classes up to row/column permutations and switchings via pynauty; phases by CP-SAT."""
import itertools, pickle, sys, numpy as np
import pynauty
from ortools.sat.python import cp_model
n = 9
vecs = [np.array((1,) + v) for v in itertools.product((1, -1), repeat=n - 1) if sum(1 for t in v if t < 0) in (0, 3, 6)]
V = len(vecs); IP = np.array([[int(a @ b) for b in vecs] for a in vecs])
ok = (np.abs(IP) == 3) | (np.abs(IP) == 9)
mats = []
def dfs(rows):
    if len(rows) == n:
        M = np.array([vecs[i] for i in rows])
        G = M.T @ M
        if all(abs(G[i, j]) in (3, 9) for i in range(n) for j in range(i + 1, n)): mats.append(M)
        return
    start = rows[-1] if len(rows) > 1 else 0
    for v in range(start, V):
        if all(ok[v, r] for r in rows): dfs(rows + [v])
allones = next(i for i, v in enumerate(vecs) if (v == 1).all())
dfs([allones])
print("normalised matrices:", len(mats), flush=True)
def canon(M):
    # signed bipartite graph: vertices rows r_i, -r_i, cols c_j, -c_j; r_i ~ c_j and -r_i ~ -c_j if M=+1, r_i ~ -c_j, -r_i ~ c_j if -1;
    # plus r_i ~ -r_i and c_j ~ -c_j (marked by colour); switching = swapping r_i <-> -r_i
    N = 4 * n; adj = {v: [] for v in range(N)}
    R = lambda i, s: i + (n if s else 0); C = lambda j, s: 2 * n + j + (n if s else 0)
    for i in range(n):
        adj[R(i, 0)].append(R(i, 1))
        for j in range(n):
            if M[i, j] > 0: adj[R(i, 0)].append(C(j, 0)); adj[R(i, 1)].append(C(j, 1))
            else: adj[R(i, 0)].append(C(j, 1)); adj[R(i, 1)].append(C(j, 0))
    for j in range(n): adj[C(j, 0)].append(C(j, 1))
    g = pynauty.Graph(N, adjacency_dict=adj, vertex_coloring=[set(range(2 * n)), set(range(2 * n, 4 * n))])
    return pynauty.certificate(g)
classes = {}
for M in mats:
    classes.setdefault(canon(M), M)
# transposes are equivalent placements up to swapping sides; keep both (a placement uses M or M^T)
print("classes up to switching and permutations:", len(classes), flush=True)
def phases(M):
    m = cp_model.CpModel()
    c = [[m.NewIntVar(0, 2, '') for j in range(n)] for i in range(n)]
    for j in range(n): m.Add(c[0][j] == 0)
    for i in range(n): m.Add(c[i][0] == 0)
    for (P, Q) in ((M, c), (M.T, [list(r) for r in zip(*c)])):
        for i in range(n):
            for k in range(i + 1, n):
                cnt = [[] for _ in range(6)]
                for j in range(n):
                    base = 0 if P[i, j] * P[k, j] > 0 else 3
                    for dv in range(-2, 3):
                        b = m.NewBoolVar('')
                        m.Add(Q[i][j] - Q[k][j] == dv).OnlyEnforceIf(b); m.Add(Q[i][j] - Q[k][j] != dv).OnlyEnforceIf(b.Not())
                        cnt[(base + 2 * dv) % 6].append(b)
                nn = [sum(x) for x in cnt]
                m.Add(nn[0] - nn[3] == nn[2] - nn[5]); m.Add(nn[2] - nn[5] == nn[4] - nn[1])
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = 120; s.parameters.num_workers = 1
    r = s.Solve(m)
    return s.StatusName(r), ([[s.Value(c[i][j]) for j in range(n)] for i in range(n)] if r in (cp_model.OPTIMAL, cp_model.FEASIBLE) else None)
good = []
for key, M in classes.items():
    st, ph = phases(M)
    G = M @ M.T
    if ph is not None: good.append((M, ph))
    print("class", len(good), st, "row Gram off-diag values", sorted(set(int(G[i, j]) for i in range(n) for j in range(i + 1, n))), flush=True)
pickle.dump(good, open("bh9_classes.pkl", "wb"))
print("BH(9,6) sign classes with phases:", len(good))

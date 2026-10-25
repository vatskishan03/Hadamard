"""All 7x7 +-1 sign patterns M (up to row/column switching and permutations) with every row and every column inner
product equal to +-3 (necessary for the sign pattern of a BH(7, 6): seven sixth roots summing to 0 are one triangle
and two antipodal pairs), and for each class whether cube-root phases make M o Omega a BH(7, 6) (CP-SAT)."""
import itertools, numpy as np
from ortools.sat.python import cp_model
vecs = [np.array((1,) + v) for v in itertools.product((1, -1), repeat=6)]
rows_all = []
def dfs(rows):
    if len(rows) == 7:
        M = np.array(rows)
        G = M.T @ M
        if all(abs(G[i, j]) == 3 for i in range(7) for j in range(i + 1, 7)): rows_all.append(M)
        return
    for v in vecs:
        if len(rows) >= 2 and tuple(v) <= tuple(rows[-1]): continue
        if all(abs(int(v @ r)) == 3 for r in rows):
            dfs(rows + [v])
# first row all ones (row switching makes first column +; column switching makes first row +)
dfs([np.ones(7, int)])
print("matrices with first row/col normalised:", len(rows_all))
# canonical form under switching + permutations: use invariant (sorted spectrum of M M^T as Seidel) and brute canon
def canon(M):
    best = None
    for perm in itertools.permutations(range(7)):
        A = M[:, perm]
        # switch rows so first column +, then sort rows; switch columns so first row + ... approximate canon via sorted Gram
        G = (A.T @ A)
        key = tuple(sorted(tuple(sorted(G[i])) for i in range(7)))
        if best is None or key < best: best = key
    return best
classes = {}
for M in rows_all:
    key = (tuple(sorted(np.linalg.eigvalsh(M @ M.T).round(4))), tuple(sorted(np.linalg.eigvalsh(M.T @ M).round(4))))
    classes.setdefault(key, []).append(M)
print("spectral classes:", len(classes))
def phases(M):
    m = cp_model.CpModel()
    c = [[m.NewIntVar(0, 2, '') for j in range(7)] for i in range(7)]
    for j in range(7): m.Add(c[0][j] == 0)
    for i in range(7): m.Add(c[i][0] == 0)
    # entry = M_ij * w^{c_ij} = zeta6^{e_ij}, e = 2c (+3 if M = -1); row pair sum zero <=> counts n0..n5 satisfy
    # n0 - n3 = n2 - n5 = n4 - n1
    for i in range(7):
        for k in range(i + 1, 7):
            cnt = [[] for _ in range(6)]
            for j in range(7):
                base = 0 if M[i, j] * M[k, j] > 0 else 3
                diff = m.NewIntVar(-2, 2, ''); m.Add(diff == c[i][j] - c[k][j])
                for dv in range(-2, 3):
                    b = m.NewBoolVar(''); m.Add(diff == dv).OnlyEnforceIf(b); m.Add(diff != dv).OnlyEnforceIf(b.Not())
                    cnt[(base + 2 * dv) % 6].append(b)
            n = [sum(x) for x in cnt]
            m.Add(n[0] - n[3] == n[2] - n[5]); m.Add(n[2] - n[5] == n[4] - n[1])
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = 60
    return s.StatusName(s.Solve(m))
for key, Ms in classes.items():
    print("Gram spectra", key, "count", len(Ms), "phases:", phases(Ms[0]), flush=True)

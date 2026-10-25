"""Find any symmetric Hadamard matrix K of order N = 4w with a symmetric 4-split K = W0 + W1 + W2 + W3 (disjoint
symmetric weighing matrices W(N, w)). Variables: label a and sign of each entry (i <= j). Conditions: each W_a^2 = wI and
K^2 = N I. usage: symsplit_any.py N seconds workers"""
import sys, time, numpy as np
from ortools.sat.python import cp_model
sys.path.insert(0, '.')
from split4 import plug
N, T, W = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]); w = N // 4
m = cp_model.CpModel()
X = {}   # X[i,j,a,s]: entry (i,j) in layer a with sign s (0: +, 1: -)
for i in range(N):
    for j in range(i, N):
        for a in range(4):
            for s in range(2): X[i, j, a, s] = m.NewBoolVar('')
        m.AddExactlyOne([X[i, j, a, s] for a in range(4) for s in range(2)])
x = lambda i, j, a, s: X[min(i, j), max(i, j), a, s]
lay = lambda i, j, a: [x(i, j, a, 0), x(i, j, a, 1)]
def prod(u, v):
    y = m.NewBoolVar(''); m.AddBoolAnd([u, v]).OnlyEnforceIf(y); m.AddBoolOr([u.Not(), v.Not()]).OnlyEnforceIf(y.Not()); return y
for a in range(4):
    for i in range(N): m.Add(sum(x(i, j, a, s) for j in range(N) for s in range(2)) == w)
for i in range(N):
    for i2 in range(i + 1, N):
        allpos, allneg = [], []
        for a in range(4):
            pos, neg = [], []
            for j in range(N):
                for s in range(2):
                    for s2 in range(2):
                        y = prod(x(i, j, a, s), x(i2, j, a, s2))
                        (pos if s == s2 else neg).append(y)
            m.Add(sum(pos) == sum(neg))
        # K rows orthogonal: entries signs regardless of layer
        for j in range(N):
            sp = [m.NewBoolVar(''), m.NewBoolVar('')]
            for s in range(2): m.Add(sum(x(i, j, a, s) for a in range(4)) == 1).OnlyEnforceIf(sp[s]); m.Add(sum(x(i, j, a, s) for a in range(4)) == 0).OnlyEnforceIf(sp[s].Not())
            sq = [m.NewBoolVar(''), m.NewBoolVar('')]
            for s in range(2): m.Add(sum(x(i2, j, a, s) for a in range(4)) == 1).OnlyEnforceIf(sq[s]); m.Add(sum(x(i2, j, a, s) for a in range(4)) == 0).OnlyEnforceIf(sq[s].Not())
            for s in range(2):
                for s2 in range(2):
                    (allpos if s == s2 else allneg).append(prod(sp[s], sq[s2]))
        m.Add(sum(allpos) == sum(allneg))
# normalise: first row all +, diagonal entry (0,0) in layer 0
for j in range(N): m.Add(sum(x(0, j, a, 0) for a in range(4)) == 1)
m.Add(x(0, 0, 0, 0) == 1)
s_ = cp_model.CpSolver(); s_.parameters.max_time_in_seconds = T; s_.parameters.num_workers = W
t0 = time.time(); r = s_.Solve(m); print(f"N={N} {s_.StatusName(r)} {time.time()-t0:.1f}s", flush=True)
if r in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    parts = [np.array([[sum((1 - 2 * s) * s_.Value(x(i, j, a, s)) for s in range(2)) for j in range(N)] for i in range(N)], np.int64) for a in range(4)]
    K = sum(parts)
    print("K sym Hadamard:", (K == K.T).all() and np.array_equal(K @ K, N * np.eye(N, dtype=np.int64)),
          "parts:", [bool((P == P.T).all() and np.array_equal(P @ P, w * np.eye(N, dtype=np.int64))) for P in parts])
    H = plug(parts); n = H.shape[0]
    print("H symmetric Hadamard of order", n, ":", bool((H == H.T).all()) and np.array_equal(H @ H, n * np.eye(n, dtype=np.int64)))
    np.save(f"symsplit_{N}.npy", np.array(parts))

"""Does a symmetric weighing matrix W(n, k) exist (W = W^T, entries 0, +-1, W^2 = k I)? CP-SAT, small n."""
import sys, time
from ortools.sat.python import cp_model
n, k, T = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
m = cp_model.CpModel()
P = {}; N = {}
for i in range(n):
    for j in range(i, n):
        P[i, j] = m.NewBoolVar(''); N[i, j] = m.NewBoolVar(''); m.AddBoolOr([P[i, j].Not(), N[i, j].Not()])
p = lambda i, j: P[min(i, j), max(i, j)]; q = lambda i, j: N[min(i, j), max(i, j)]
for i in range(n): m.Add(sum(p(i, j) + q(i, j) for j in range(n)) == k)
for i in range(n):
    for i2 in range(i + 1, n):
        pos, neg = [], []
        for j in range(n):
            for (a, b, sgn) in ((p(i, j), p(i2, j), 1), (q(i, j), q(i2, j), 1), (p(i, j), q(i2, j), -1), (q(i, j), p(i2, j), -1)):
                y = m.NewBoolVar(''); m.AddBoolAnd([a, b]).OnlyEnforceIf(y); m.AddBoolOr([a.Not(), b.Not()]).OnlyEnforceIf(y.Not())
                (pos if sgn > 0 else neg).append(y)
        m.Add(sum(pos) == sum(neg))
s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = T; s.parameters.num_workers = 1
r = s.Solve(m); print(f"symmetric W({n},{k}): {s.StatusName(r)}", flush=True)
if r in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    for i in range(n): print(' '.join('+' if s.Value(p(i, j)) else '-' if s.Value(q(i, j)) else '0' for j in range(n)))

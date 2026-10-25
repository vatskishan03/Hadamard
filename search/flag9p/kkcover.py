"""Is there a decomposition of K_p into four cliques K_{d+1} and six bicliques K_{d,d} (d = (p-1)/4) with the D-balance of
the flag equation (NOTES.md §96): biclique sides have D = +1 and -1, a clique has (d+1)/2 flags of each sign, a point
has 2 flags of each sign? CP-SAT. usage: kkcover.py p seconds"""
import sys, time
from ortools.sat.python import cp_model
p, T = int(sys.argv[1]), float(sys.argv[2]); d = (p - 1) // 4
m = cp_model.CpModel()
q = [[m.NewBoolVar('') for x in range(p)] for c in range(4)]
dq = [[m.NewBoolVar('') for x in range(p)] for c in range(4)]      # D = +1 for x's flag in clique c (if member)
L = [[m.NewBoolVar('') for x in range(p)] for k in range(6)]
R = [[m.NewBoolVar('') for x in range(p)] for k in range(6)]
def AND(a, b):
    y = m.NewBoolVar(''); m.AddBoolAnd([a, b]).OnlyEnforceIf(y); m.AddBoolOr([a.Not(), b.Not()]).OnlyEnforceIf(y.Not()); return y
for c in range(4): m.Add(sum(q[c]) == d + 1)
for k in range(6):
    m.Add(sum(L[k]) == d); m.Add(sum(R[k]) == d)
    for x in range(p): m.AddBoolOr([L[k][x].Not(), R[k][x].Not()])
for x in range(p):
    for y in range(x + 1, p):
        terms = [AND(q[c][x], q[c][y]) for c in range(4)]
        terms += [AND(L[k][x], R[k][y]) for k in range(6)] + [AND(R[k][x], L[k][y]) for k in range(6)]
        m.Add(sum(terms) == 1)
# D-balance
pos = [[AND(q[c][x], dq[c][x]) for x in range(p)] for c in range(4)]
for c in range(4): m.Add(sum(pos[c]) == (d + 1) // 2)
for x in range(p):
    m.Add(sum(pos[c][x] for c in range(4)) + sum(L[k][x] for k in range(6)) == 2)
    m.Add(sum(q[c][x] for c in range(4)) + sum(L[k][x] + R[k][x] for k in range(6)) == 4)
# symmetry: point 0 in clique 0
m.Add(q[0][0] == 1)
s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = T; s.parameters.num_workers = 2
t0 = time.time(); r = s.Solve(m); print(f"p={p} d={d}: {s.StatusName(r)} ({time.time() - t0:.0f}s)", flush=True)
if r in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    for c in range(4): print("clique", c, [x for x in range(p) if s.Value(q[c][x])], "D+", [x for x in range(p) if s.Value(q[c][x]) and s.Value(dq[c][x])])
    for k in range(6): print("biclique", k, [x for x in range(p) if s.Value(L[k][x])], [x for x in range(p) if s.Value(R[k][x])])

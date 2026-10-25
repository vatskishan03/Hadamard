"""H-invariant circulant 4-split of a Goethals-Seidel quadruple of order 191 (dk764.py) for the x9 plug-in:
labels on (sequence, slot) (slot = {0} or an H-coset), each layer: sum_k PAF(layer_k) = 191 delta.
usage: split1719.py which seconds workers seed"""
import sys, time, json, numpy as np
from ortools.sat.python import cp_model
sys.path.insert(0, '.')
from dk764 import n, quad, slots
which, T, W, seed = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
Q = quad(which); sl = slots(); NS = 39
size = [1] + [5] * 38
val = [[int(Q[k][np.where(sl == s)[0][0]]) for s in range(NS)] for k in range(4)]
# orbit representatives of nonzero shifts under <H, -1>
seen = set(); dreps = []
for d in range(1, n):
    if d in seen: continue
    dreps.append(d)
    for h in [1, 39, 184, 109, 49]:
        seen.add(d * h % n); seen.add((-d * h) % n)
assert len(dreps) == 19
# coefficient of [both slots s1, s2 in layer a] in PAF_k(d): sum over t of x(t) x(t+d)
coef = {}
for k in range(4):
    for d in dreps:
        c = {}
        for t in range(n):
            s1, s2 = int(sl[t]), int(sl[(t + d) % n])
            key = (min(s1, s2), max(s1, s2))
            c[key] = c.get(key, 0) + int(Q[k][t] * Q[k][(t + d) % n])
        coef[k, d] = {key: v for key, v in c.items() if v}
m = cp_model.CpModel()
L = {(k, s, a): m.NewBoolVar(f'L{k}_{s}_{a}') for k in range(4) for s in range(NS) for a in range(4)}
for k in range(4):
    for s in range(NS): m.AddExactlyOne([L[k, s, a] for a in range(4)])
Y = {}
def y(k, s1, s2, a):
    if s1 == s2: return L[k, s1, a]
    key = (k, s1, s2, a)
    if key not in Y:
        v = m.NewBoolVar(''); Y[key] = v
        m.AddBoolAnd([L[k, s1, a], L[k, s2, a]]).OnlyEnforceIf(v)
        m.AddBoolOr([L[k, s1, a].Not(), L[k, s2, a].Not()]).OnlyEnforceIf(v.Not())
    return Y[key]
for a in range(4):
    m.Add(sum(size[s] * L[k, s, a] for k in range(4) for s in range(NS)) == n)
    for d in dreps:
        m.Add(sum(v * y(k, s1, s2, a) for k in range(4) for (s1, s2), v in coef[k, d].items()) == 0)
# label symmetry: the zero slot of sequence k is in layer k (each layer holds exactly one zero slot)
for k in range(4): m.Add(L[k, 0, k] == 1)
s_ = cp_model.CpSolver(); s_.parameters.max_time_in_seconds = T; s_.parameters.num_workers = W
s_.parameters.random_seed = seed; s_.parameters.log_search_progress = False
t0 = time.time(); r = s_.Solve(m)
print(f"which={which} seed={seed} status={s_.StatusName(r)} time={time.time() - t0:.1f}s aux={len(Y)}", flush=True)
if r in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    lab = [[next(a for a in range(4) if s_.Value(L[k, s, a])) for s in range(NS)] for k in range(4)]
    json.dump({"which": which, "labels": lab}, open(f"split_{which}_{seed}.json", "w"))
    print("labels", lab, flush=True)

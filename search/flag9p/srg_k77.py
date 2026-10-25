"""Does a conference graph SRG(4d+1, 2d, d-1, d) contain disjoint d-sets A, B whose Seidel signs S[A, B] (S = J - I - 2A)
equal u_a v_b (1 - 2 [b = pi(a)]) (the sign class J - 2P of BH(d, 6))? A = {0..d-1}, B = {d..2d-1}, pi = identity (by
relabelling). CP-SAT. usage: srg_k77.py p seconds workers"""
import sys, time
from ortools.sat.python import cp_model
p, T, W = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]); d = (p - 1) // 4
k, lam, mu = 2 * d, d - 1, d
m = cp_model.CpModel()
a = {}
for x in range(p):
    for y in range(x + 1, p): a[x, y] = m.NewBoolVar('')
E = lambda x, y: a[min(x, y), max(x, y)]
for x in range(p): m.Add(sum(E(x, y) for y in range(p) if y != x) == k)
for x in range(p):
    for y in range(x + 1, p):
        prods = []
        for z in range(p):
            if z in (x, y): continue
            t = m.NewBoolVar(''); m.AddBoolAnd([E(x, z), E(z, y)]).OnlyEnforceIf(t)
            m.AddBoolOr([E(x, z).Not(), E(z, y).Not()]).OnlyEnforceIf(t.Not()); prods.append(t)
        m.Add(sum(prods) == lam).OnlyEnforceIf(E(x, y)); m.Add(sum(prods) == mu).OnlyEnforceIf(E(x, y).Not())
U = [m.NewBoolVar('') for _ in range(d)]; V = [m.NewBoolVar('') for _ in range(d)]
m.Add(U[0] == 0)
for i in range(d):
    for j in range(d):
        # S = -1 (edge)  <=>  U_i xor V_j xor [i = j] = 1
        if i == j: m.AddBoolXOr([E(i, d + j), U[i], V[j]])            # edge = not(U xor V)
        else: m.AddBoolXOr([E(i, d + j), U[i], V[j], m.NewConstant(1)])
s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = T; s.parameters.num_workers = W
t0 = time.time(); r = s.Solve(m); print(f"p={p}: SRG with a J-2P K_(d,d): {s.StatusName(r)} ({time.time() - t0:.0f}s)", flush=True)
if r in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    adj = [[1 if x != y and s.Value(E(x, y)) else 0 for y in range(p)] for x in range(p)]
    open(f"srg{p}_k77.txt", "w").write('\n'.join(''.join(map(str, r_)) for r_ in adj) + '\n')
    print("saved", flush=True)

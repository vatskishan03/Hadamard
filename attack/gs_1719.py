import itertools, sys, time
import numpy as np
sys.path.insert(0, '.')
import spence_lib as L, gs_lib as GS

def bd(*b):
    n = sum(x.shape[0] for x in b); M = np.zeros((n, n), dtype=np.int64); i = 0
    for x in b: k = x.shape[0]; M[i:i+k, i:i+k] = x; i += k
    return M
I2 = np.eye(2, dtype=np.int64); I1 = np.eye(1, dtype=np.int64)
G = L.AbGroup([3, 3, 191])
# automorphisms of Z_3^2
U  = np.array([[1, 1], [0, 1]])      # unipotent, order 3
D  = np.array([[1, 0], [0, 2]])      # reflection diag(1,-1), order 2
D2 = np.array([[2, 0], [0, 1]])      # diag(-1, 1)
S9 = L.field_mult_matrix(3, 2)       # Singer, order 8
mats2 = {"1": I2, "U": U, "D": D, "D2": D2, "UD": None, "S^2": L.mat_pow_mod(S9, 2, 3), "S": S9}
# Z_191 multipliers of each order
def gen_of_order(p, d):
    for g in range(2, p):
        if pow(g, d, p) == 1 and all(pow(g, d // q, p) != 1 for q in {f for f in range(2, d + 1) if d % f == 0 and all(f % r for r in range(2, f))}):
            return g
units = {d: gen_of_order(191, d) for d in (1, 2, 5, 10, 19, 38, 95, 190)}
units[1] = 1
rows = []
two_part = [("1", [I2]), ("U", [U]), ("D", [D]), ("B=<U,D>", [U, D]), ("<U,D2>", [U, D2]), ("S^2", [L.mat_pow_mod(S9,2,3)]), ("S", [S9])]
for (n2, g2s), d in itertools.product(two_part, (1, 2, 5, 10, 19, 38, 95, 190)):
    gens = [bd(M, I1) for M in g2s] + ([bd(I2, units[d] * I1)] if d > 1 else [])
    if not gens: continue
    st = L.Structure(G, [G.perm_from_matrix(M) for M in gens])
    if not (8 <= st.O <= 64): continue
    t = time.time(); b = GS.gs_estimate(st, samples=2500)
    if b is None:
        rows.append((99, n2, d, st.O, st.C, None)); continue
    rows.append((b['lW'], n2, d, st.O, st.C, b))
for lW, n2, d, O, C, b in sorted(rows, key=lambda r: r[0]):
    if b is None: print(f"  Z3^2:{n2:8s} Z191:ord {d:3d}  O={O:3d} C={C:3d}  no admissible pattern"); continue
    print(f"  Z3^2:{n2:8s} Z191:ord {d:3d}  O={O:3d} C={C:3d}  pattern={b['pattern']}  logE={b['lE']:6.1f} logN={b['lN']:5.1f} logW(L=1e8)={b['lW']:5.1f}")

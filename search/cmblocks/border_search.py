# Search for 4x4 border data (A, B, B') for the quaternary Miyamoto array such that
#  H0: A A* = 4I, B B* = 4I, B' B'* = 4I
#  H1: A B'* + 2 B R* = 0                (rows of the array orthogonal)
#  C1: A - 2R = lam * B = lam' * B'      (all sixteen bordered blocks commute)
#  T1: B' = kappa * B                    (one (anti-linear) monomial R transposes all blocks)
# R = block row sums of the Hermitian core (0, i / -i, 0 pattern), possibly with the variants of the core.
import itertools, numpy as np, sys
U4 = [1, 1j, -1, -1j]
def variants():
    out = {}
    for e1 in (1, -1):
        for e2 in (1, -1):
            R = np.zeros((4, 4), complex)
            R[0, 1], R[1, 0] = 1j * e1, -1j * e1
            R[2, 3], R[3, 2] = 1j * e2, -1j * e2
            out[(e1, e2)] = R
    return out
def search(R, verbose=False):
    sols = []
    nz = [(r, c) for r in range(4) for c in range(4) if abs(R[r, c]) > 0]
    free = [(r, c) for r in range(4) for c in range(4) if (r, c) not in nz]
    for lam in U4:
        for kappa in U4:
            # B on nz entries: lam B = -R
            Bfix = {rc: -R[rc] / lam for rc in nz}
            # enumerate rows of B incrementally with row-orthogonality BB* = 4I
            rows = []
            def rowopts(r):
                opts = []
                cols = [c for c in range(4)]
                fr = [c for c in cols if (r, c) in free]
                for vals in itertools.product(U4, repeat=len(fr)):
                    row = np.zeros(4, complex)
                    for c in cols:
                        if (r, c) in Bfix: row[c] = Bfix[(r, c)]
                    for c, v in zip(fr, vals): row[c] = v
                    opts.append(row)
                return opts
            opts = [rowopts(r) for r in range(4)]
            def rec(r, chosen):
                if r == 4:
                    B = np.array(chosen)
                    A = 2 * R + lam * B
                    if not np.allclose(np.abs(A), 1): return
                    if not np.allclose(A @ A.conj().T, 4 * np.eye(4)): return
                    Bp = kappa * B
                    if not np.allclose(A @ Bp.conj().T + 2 * B @ R.conj().T, 0): return
                    sols.append((lam, kappa, B.copy(), A.copy()))
                    return
                for row in opts[r]:
                    if all(abs(np.vdot(prev, row)) < 1e-9 for prev in chosen):
                        rec(r + 1, chosen + [row])
            rec(0, [])
    return sols
if __name__ == "__main__":
    for key, R in variants().items():
        s = search(R)
        print("R variant", key, "solutions:", len(s))
        for lam, kappa, B, A in s[:3]:
            print(" lam", lam, "kappa", kappa); print(np.round(B, 2)); print(np.round(A, 2))

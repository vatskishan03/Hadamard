# Complex (linear) Goethals-Seidel plug-in on the quaternary blocks: which identity fails?
import numpy as np, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import cmb, paley_circ as P
Kr, Ki, m = P.block_circulant_K(11)
q = 2 * m + 1
N = cmb.build_blocks(q, cmb.prim_root(q), Kr, Ki)
Rb = cmb.reversal_perm(q)
Rm = np.eye(q)[Rb].T if False else np.zeros((q, q)); Rm[np.arange(q), Rb] = 1
T = cmb.t_matrices(3); t = 3
Ct = [cmb.circ_int(r) for r in T]
A = [sum(np.kron(Ct[k], N[k][j]) for k in range(4)) for j in range(4)]
Rt = np.zeros((t, t)); Rt[np.arange(t), (-np.arange(t)) % t] = 1
R = np.kron(Rt, Rm)
print("complex A_j commute:", all(np.allclose(x @ y, y @ x) for x in A for y in A),
      "  R A R = A^T:", all(np.allclose(R @ x @ R, x.T) for x in A),
      "  sum A_j A_j* = 4qt I:", np.allclose(sum(x @ x.conj().T for x in A), 4 * q * t * np.eye(q * t)))
B_, C_, D_ = A[1], A[2], A[3]
AR = [x @ R for x in A]; ATR = [x.T @ R for x in A]
Z = lambda: None
G = np.block([[A[0], AR[1], AR[2], AR[3]], [-AR[1], A[0], ATR[3], -ATR[2]], [-AR[2], -ATR[3], A[0], ATR[1]], [-AR[3], ATR[2], -ATR[1], A[0]]])
n = G.shape[0]
M = G @ G.conj().T
off = np.abs(M - 4 * q * t * np.eye(n)).max()
print(f"complex GS array (q={q}, t=3): order {n}, entries units: {np.allclose(np.abs(G), 1)}, max |G G* - 4qt I| = {off:.1f}")
n1 = q * t
blk12 = M[:n1, n1:2 * n1]
pred = (-A[0] @ R @ B_.conj().T + B_ @ R @ A[0].conj().T) + (C_ @ C_.conj().T * 0)  # placeholder
# the cross term predicted by the algebra: (-A conj(B) + B conj(A)) R + C conj(D) - D conj(C)
cA, cB, cC, cD = [x.conj() for x in A]
pred = (-A[0] @ cB + B_ @ cA) @ R + C_ @ cD - D_ @ cC
print("block (1,2) of G G* equals (-A conj(B) + B conj(A)) R + C conj(D) - D conj(C):", np.allclose(blk12, pred),
      "  its norm:", np.abs(blk12).max())
print("  A conj(B) = B conj(A)?", np.allclose(A[0] @ cB, B_ @ cA), "  C conj(D) = D conj(C)?", np.allclose(C_ @ cD, D_ @ cC))

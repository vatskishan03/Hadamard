# Do the sixteen blocks of the existing complex_miyamoto (borders B' = R P F, A = -B P F/2) commute / satisfy R N R = N^T ?
import numpy as np, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import cmb, paley_circ as P
Kr, Ki, m = P.block_circulant_K(11)
q = 2 * m + 1
g = cmb.prim_root(q)
U, V = cmb.core_blocks(q, g, Kr, Ki)
R = np.array([[0, 1j, 0, 0], [-1j, 0, 0, 0], [0, 0, 0, 1j], [0, 0, -1j, 0]])
B = np.array([[1, 1, 1, 1], [1, 1, -1, -1], [1, -1, 1, -1], [1, -1, -1, 1]], complex)
F = np.array([[1, 1, 1, 1], [1, 1, -1, -1], [1, -1, 1, -1], [-1, 1, 1, -1]], complex)
Pm = np.eye(4)[[1, 0, 3, 2]]
Bp, A = R @ Pm @ F, -B @ Pm @ F / 2
N = [[None] * 4 for _ in range(4)]
for r in range(4):
    for c in range(4):
        Z = np.zeros((q, q), complex)
        Z[0, 0] = A[r, c]; Z[0, 1:] = B[r, c]; Z[1:, 0] = Bp[r, c]
        u, v = U[r][c], V[r][c]
        Z[1:m + 1, 1:m + 1] = u + v; Z[1 + m:, 1 + m:] = u + v; Z[1:m + 1, 1 + m:] = u - v; Z[1 + m:, 1:m + 1] = u - v
        N[r][c] = Z
H = np.block(N)
print("old borders: H H* = 4qI:", np.allclose(H @ H.conj().T, 4 * q * np.eye(4 * q)))
flat = [N[r][c] for r in range(4) for c in range(4)]
bad = sum(1 for i, x in enumerate(flat) for y in flat[i + 1:] if not np.allclose(x @ y, y @ x))
print("old borders: non-commuting pairs among the 16 blocks:", bad, "of", 120)
Rb = cmb.reversal_perm(q)
print("old borders: blocks with R N R = N^T:", sum(np.array_equal(x[np.ix_(Rb, Rb)], x.T) for x in flat), "of 16")
print("old borders: B'/B entrywise (should be constant for commuting blocks):"); print(np.round(Bp / B, 2))
print("old borders: (A - 2R)/B entrywise:"); print(np.round((A - 2 * R) / B, 2))

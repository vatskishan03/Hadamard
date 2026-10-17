"""Prepare data for the d=3 walk search over Z_1893: (X0, X1, mu X1, mu^2 X1), H=<8>, K=<8,2>, mu=2.
Writes:  pool0.bin  int32 N0, then int16 V0[N0][15]   (class values of every X0 with row sum 87)
         pool0_y.bin int8 Y0[N0][30]                 (orbit signs of each X0, K-orbit order)
         walk.bin   int32 k1, CK; int32 sizes[k1]; int32 MF[CK][k1][k1]; int32 DW[CK]; double mu0[CK], sd0[CK]
The class value of X1 on K-class k is W_k = sum_j PAF_X1(mu^j s_k) = DW_k + 1/2 y^T MF_k y."""
import os
import sys, struct
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import d3_model as D
n, hg, mu = 1893, 8, 2
G, stH, stK, fold = D.setup(n, hg, mu)
CK, k0, k1 = stK.C, stK.O, stH.O
# ---- X0 pool: orbit order by size -----------------------------------------------------
S0 = stK.S                                               # (CK, k0, k0)
sizes0 = stK.sizes
small = [o for o in range(k0) if sizes0[o] <= 2]
o45 = [o for o in range(k0) if sizes0[o] == 45]; o90 = [o for o in range(k0) if sizes0[o] == 90]
def signs(m):
    idx = np.arange(1 << m)
    return (1 - 2 * ((idx[:, None] >> np.arange(m)) & 1)).astype(np.int8)
Y45, Y90 = signs(len(o45)), signs(len(o90))
s45, s90 = Y45.sum(1), Y90.sum(1)
ysmall = {o: (-1 if sizes0[o] == 1 else -1) for o in small}       # y0 = -1, y_(size 2) = -1 forced by row sum 87
rows, V0s = [], []
M0 = S0 + S0.transpose(0, 2, 1)
for c in range(CK): np.fill_diagonal(M0[c], 0)
D0 = np.array([np.trace(S0[c]) for c in range(CK)])
for b in range(-14, 15, 2):
    a = 2 * (1 - b)
    if abs(a) > 14: continue
    A = Y45[s45 == a]; B = Y90[s90 == b]
    if len(A) == 0 or len(B) == 0: continue
    # full y = small + A part + B part ; V0_c = D0_c + 1/2 y^T M0_c y
    for i0 in range(0, len(A), 2048):
        Ai = A[i0:i0 + 2048]
        nA, nB = len(Ai), len(B)
        Y = np.zeros((nA * nB, k0), dtype=np.int8)
        for o, v in ysmall.items(): Y[:, o] = v
        Y[:, o45] = np.repeat(Ai, nB, axis=0); Y[:, o90] = np.tile(B, (nA, 1))
        Yf = Y.astype(np.float32)
        V = D0[None, :] + 0.5 * np.einsum('bo,cop,bp->bc', Yf, M0.astype(np.float32), Yf, optimize=True)
        rows.append(Y); V0s.append(np.rint(V).astype(np.int16))
Y0 = np.concatenate(rows); V0 = np.concatenate(V0s)
assert np.all((Y0.astype(int) @ sizes0) == 87)
print('X0 pool', len(Y0), 'row sums ok')
# spot check against direct PAF
x0 = Y0[123][stK.orb].astype(int)
for k, s in enumerate(stK.reps[:3]): assert int(x0 @ np.roll(x0, -s)) == int(V0[123, k]), 'V0 mismatch'
with open('pool0.bin', 'wb') as f: f.write(struct.pack('i', len(V0))); f.write(V0.tobytes())
Y0.tofile('pool0_y.bin')
# ---- X1 walk data ------------------------------------------------------------------------
S1 = stH.S; M1 = S1 + S1.transpose(0, 2, 1)
for c in range(stH.C): np.fill_diagonal(M1[c], 0)
D1 = np.array([np.trace(S1[c]) for c in range(stH.C)])
MF = np.stack([M1[fold[k]].sum(0) for k in range(CK)]).astype(np.int32)
DW = np.array([D1[fold[k]].sum() for k in range(CK)], dtype=np.int32)
mu0 = V0.astype(np.float64).mean(0); sd0 = V0.astype(np.float64).std(0)
with open('walk.bin', 'wb') as f:
    f.write(struct.pack('ii', k1, CK)); f.write(stH.sizes.astype(np.int32).tobytes())
    f.write(MF.tobytes()); f.write(DW.tobytes()); f.write(mu0.tobytes()); f.write(sd0.tobytes())
np.save('orb.npy', np.stack([stK.orb, stH.orb]))
print('walk data: k1', k1, 'CK', CK, 'mu0', mu0.round(1).tolist(), 'sd0', sd0.round(1).tolist())

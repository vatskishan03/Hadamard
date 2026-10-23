"""independent exact check (own numpy code, no repo code) of a 'goethals-seidel-group' skew certificate over Z_m x Z_q (mixed radix k*q + f)"""
import json, sys, numpy as np
c = json.load(open(sys.argv[1]))
mods = c['group']; n = int(np.prod(mods)); assert len(mods) == 2
q = mods[1]
seqs = [np.array([1 if ch == '+' else -1 for ch in s], dtype=np.int64) for s in c['sequences']]
def idx(k, f): return (k % mods[0]) * q + (f % q)
def dev(x):
    X = np.zeros((n, n), dtype=np.int64)
    for g in range(n):
        kg, fg = divmod(g, q)
        for h in range(n):
            kh, fh = divmod(h, q); X[g, h] = x[idx(kh - kg, fh - fg)]
    return X
R = np.zeros((n, n), dtype=np.int64)
for g in range(n):
    kg, fg = divmod(g, q); R[g, idx(-kg, -fg)] = 1
A, B, C, D = [dev(x) for x in seqs]
H = np.block([[A, B @ R, C @ R, D @ R], [-B @ R, A, D.T @ R, -C.T @ R], [-C @ R, -D.T @ R, A, B.T @ R], [-D @ R, C.T @ R, -B.T @ R, A]])
N = 4 * n; I = np.eye(N, dtype=np.int64)
a = seqs[0]
sq = sorted({x * x % q for x in range(1, q)})
res = dict(order=N, pm1=bool((np.abs(H) == 1).all()), HHT_nI=bool(((H.astype(np.float64) @ H.T.astype(np.float64) == N * I)).all()), skew_H_plus_HT_2I=bool((H + H.T == 2 * I).all()),
           A0=int(a[0]), A_skew=all(a[idx(-divmod(g, q)[0], -divmod(g, q)[1])] == -a[g] for g in range(1, n)),
           H_invariant=all(all(x[idx(k, s * f)] == x[idx(k, f)] for k in range(mods[0]) for f in range(q) for s in sq) for x in seqs), rowsums=[int(x.sum()) for x in seqs])
print(json.dumps(res))

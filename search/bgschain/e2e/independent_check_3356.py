"""Independent re-assembly of the order-3356 matrix (own code, not the repository's builder) from a certificate (or SOLUTION) JSON that holds the four
sequences; the border is recomputed by z7search for the pattern read off the sequences.  Exact Gram check in float32 (entries and sums are integers below 2^24).
    python3 e2e/independent_check_3356.py [results/certificate_3356.json]"""
import json, os, sys, hashlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import border8
import cert_lib as CL
V = CL.V
path = sys.argv[1] if len(sys.argv) > 1 else os.path.join('results', 'certificate_3356.json')
sol = json.load(open(path))
mods = sol['group']; v = 832; k = 8
seqs = [np.array([1 if c == '+' else -1 for c in s], dtype=np.int64) for s in sol['sequences']]
coords = np.array(np.unravel_index(np.arange(v), mods)).T
idx = lambda c: np.ravel_multi_index(tuple(np.moveaxis(c % mods, -1, 0)), mods)
cos = (coords[:, 0] % 2) * 4 + (coords[:, 1] % 2) * 2 + (coords[:, 2] % 2)          # cosets of S = 2Z4^3 x Z13 in smallest-element order
ys = [[int(x[cos == c].sum()) for c in range(k)] for x in seqs]
a = [y[1] for y in ys]; d = [y[0] - y[1] for y in ys]
K, E, F = border8.make_border(a, d)
one = idx(coords[None, :, :] - coords[:, None, :])        # h - g
two = idx(coords[:, None, :] + coords[None, :, :])        # g + h
neg = idx(-coords)
A = seqs[0][one]; B, C, D = (x[two] for x in seqs[1:]); Br, Cr, Dr = (x[neg][two] for x in seqs[1:])
M = np.block([[A, B, C, D], [-B, A, Dr, -Cr], [-C, -Dr, A, Br], [-D, Cr, -Br, A]])
Ex = np.concatenate([E[:, i * k + cos] for i in range(4)], axis=1)                   # 28 x 4v
Fx = np.concatenate([F[i * k + cos, :] for i in range(4)], axis=0)                   # 4v x 28
H = np.block([[K, Ex], [Fx, M]]).astype(np.float32)
n = H.shape[0]
assert n == 3356 and np.all(np.abs(H) == 1)
G = H @ H.T
print('independent assembly: order', n, ' H H^T == 3356 I :', bool(np.array_equal(G, n * np.eye(n, dtype=np.float32))))
sha = hashlib.sha256()
for ch in V.csv_chunks(H.astype(np.int8)): sha.update(ch)
print('sha256 of the canonical CSV:', sha.hexdigest(), ' (file:', sol.get('sha256'), ')')

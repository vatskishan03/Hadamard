"""Verify an x9 split of Đoković's quadruple: layers -> GS arrays (W(764, 191), disjoint, sum = K), then
H = sum_a K_a (x) Y_a of order 6876 and H H^T = 6876 I (exact: float64 products of +-1 entries, |sums| <= 6876 < 2^53)."""
import sys, json, time, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__file__) or '.')
from dk764 import n, quad, slots
from split4 import Y, DIRS
def circ(x): return np.array([np.roll(x, k) for k in range(len(x))]).astype(np.int64)
def rev(X): return X[:, ::-1]
def gs_array(A, B, C, D):
    A, B, C, D = map(circ, (A, B, C, D))
    return np.block([[A, rev(B), rev(C), rev(D)], [-rev(B), A, rev(D.T), -rev(C.T)],
                     [-rev(C), -rev(D.T), A, rev(B.T)], [-rev(D), rev(C.T), -rev(B.T), A]])
d = json.load(open(sys.argv[1])); which, lab = d["which"], d["labels"]
Q = quad(which); sl = slots()
layers = [[np.where(np.array(lab[k])[sl] == a, Q[k], 0) for k in range(4)] for a in range(4)]
for a in range(4):
    F = sum(np.abs(np.fft.fft(x)) ** 2 for x in layers[a]); paf = np.rint(np.fft.ifft(F).real).astype(np.int64)
    print("layer", a, "weight", sum(int(np.abs(x).sum()) for x in layers[a]), "paf0", paf[0], "max|paf|", int(np.abs(paf[1:]).max()))
K = gs_array(*Q); parts = [gs_array(*layers[a]) for a in range(4)]
I764 = np.eye(4 * n, dtype=np.int64)
print("K Hadamard:", np.array_equal(K @ K.T, 4 * n * I764), "sum of parts = K:", np.array_equal(sum(parts), K),
      "disjoint:", int(sum(np.abs(P) for P in parts).max()) == 1,
      "parts W(764,191):", [np.array_equal(P @ P.T, n * I764) for P in parts])
t0 = time.time()
Hm = sum(np.kron(P, Y(DIRS[a])) for a, P in enumerate(parts)).astype(np.float64)
G = Hm @ Hm.T
N = Hm.shape[0]
ok = bool(np.array_equal(np.abs(Hm), np.ones_like(Hm))) and bool(np.array_equal(G, N * np.eye(N)))
print(f"H order {N}: entries +-1 and H H^T = {N} I: {ok} ({time.time() - t0:.0f} s)")
if ok and len(sys.argv) > 2:
    np.save(sys.argv[2], Hm.astype(np.int8))

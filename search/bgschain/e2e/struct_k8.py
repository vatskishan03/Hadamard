"""usage: python3 e2e/struct_k8.py [P1|P2|P3]   (writes results/struct_k8_<pattern>.json under the current directory)
STRUCTURAL k = 8 non-split control (G = Z4^3 x Z13, order 3356) with the repository's own assembly code (verify.build_bordered_gs, with its autocorrelation
check disabled because the sequences are random, NOT a solution): for random +-1 sequences with the prescribed coset sums (pattern a, d) the border blocks of
H H^T must be exactly those of 3356 I:  K K^T + R R^T = 3356 I and K C^T + R M^T = 0  (they depend only on the coset sums).  Also checks the (2,2) block formula."""
import sys, os, json, inspect
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import cert_lib as CL, border8
V = CL.V
mods = [4, 4, 4, 13]; S_gens = [[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, 2, 0], [0, 0, 0, 1]]
h = [[2, 2, 1, 0], [1, 0, 0, 0], [2, 3, 3, 0], [0, 0, 0, 3]]
pats = {'P1': ([2, 4, 4, 4], [4, -2, -2, 2]), 'P2': ([2, 4, 4, 4], [-4, -2, 2, 2]), 'P3': ([0, 4, 4, 4], [2, -2, 2, 4])}
nm = sys.argv[1] if len(sys.argv) > 1 else 'P1'
a, d = pats[nm]; A = [a[i] + d[i] for i in range(4)]
v = 832; k = 8
coords = np.array(np.unravel_index(np.arange(v), mods)).T
cos = ((coords[:, 0] % 2) * 4 + (coords[:, 1] % 2) * 2 + (coords[:, 2] % 2))
rng = np.random.default_rng(5)
seqs = []
for i in range(4):
    x = np.zeros(v, dtype=np.int64)
    for c in range(k):
        idx = np.flatnonzero(cos == c); y = A[i] if c == 0 else a[i]
        npl = (len(idx) + y) // 2; sel = rng.permutation(len(idx)); x[idx] = -1; x[idx[sel[:npl]]] = 1
    seqs.append(x)
assert all(int(x[cos == 0].sum()) == A[i] and all(int(x[cos == c].sum()) == a[i] for c in range(1, 8)) for i, x in enumerate(seqs))
K, E, F = border8.make_border(a, d)
rec = dict(order=4 * (v + k - 1), group=mods, subgroup_generators=S_gens, gens=[h], sequences=[CL.sstr(x) for x in seqs])
cert = CL.make_certificate(rec, K, E, F, -28, 4); cert['symmetries'] = []      # random sequences are not H-invariant
src = inspect.getsource(V.build_bordered_gs).replace('raise ValueError("sequences do not satisfy the autocorrelation condition")', 'pass')
ns = dict(V.__dict__); exec(src, ns)
H, _ = ns['build_bordered_gs'](cert)
n = H.shape[0]; print('assembled', H.shape, 'entries +-1:', bool(np.all(np.abs(H) == 1)))
Hf = H.astype(np.float32)
m = 28
G11 = Hf[:m] @ Hf[:m].T; G12 = Hf[:m] @ Hf[m:].T
ok11 = np.array_equal(G11, n * np.eye(m, dtype=np.float32)); ok12 = not np.any(G12)
print(f'border blocks of H H^T (pattern {nm}): (1,1) = {n} I : {ok11};  (1,2) = 0 : {ok12}')
# (2,2) block = C C^T + M M^T ; M M^T = I4 (x) SS with SS[g,h] = sum_i PAF_i(h-g)
Cb = Hf[m:, :m]; M = Hf[m:, m:]
G22 = Cb @ Cb.T + M @ M.T
idx = lambda c: np.ravel_multi_index(tuple(np.moveaxis(c % mods, -1, 0)), mods)
paf = np.zeros(v)
for s in range(v):
    sh = idx(coords + coords[s]); paf[s] = sum(int(x @ x[sh]) for x in seqs)
SS = paf[idx(coords[None, :, :] - coords[:, None, :])]          # SS[g, h] = paf(h - g)
blk = G22[:v, :v] - (Cb @ Cb.T)[:v, :v]
okM = np.array_equal(M @ M.T, np.kron(np.eye(4, dtype=np.float32), SS.astype(np.float32)))
print('M M^T = I4 (x) (sum PAF) :', okM)
res = dict(pattern=nm, border_blocks_ok=bool(ok11 and ok12), MMt_formula_ok=bool(okM))
os.makedirs('results', exist_ok=True)
json.dump(res, open(os.path.join('results', f'struct_k8_{nm}.json'), 'w')); print(res)

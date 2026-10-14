"""Verify a lift solution: python verify573.py <wit.npy> <maskA> <maskB> <tag>"""
import sys, json, hashlib, numpy as np
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'verify'))
import verify
from lift_poly import Poly, build, make
from ss_core import paf_sum_check
S = build(3); LQ = make(S, 'Z2^2xZ13')
w = np.load(sys.argv[1]); a, b = w[:LQ.RK], w[LQ.RK:]
P = Poly(S, LQ, a, b)
mA, mB = int(sys.argv[2]), int(sys.argv[3]); tag = sys.argv[4]
sA = np.array([1 - 2 * ((mA >> i) & 1) for i in range(P.nA)]); sB = np.array([1 - 2 * ((mB >> i) & 1) for i in range(P.nB)])
ke, XA, XB = P.keys_exact(sA, sB)
print('class keys (exact):', sorted(set(ke.tolist())))
seqs = S.sequences(XA, XB)
print('entries +-1:', all(set(np.unique(x).tolist()) <= {-1, 1} for x in seqs), 'row sums', [int(x.sum()) for x in seqs])
tot = paf_sum_check(S.G, seqs)
ok = tot[0] == 4 * S.G.v and np.all(tot[1:] == -4)
print('Spence condition (direct PAF over all 572 shifts):', bool(ok), 'values', sorted(set(tot[1:].tolist())))
H4 = [[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]]
Pm = (np.ones((4, 4), int) - 2 * np.eye(4, dtype=int)).tolist()
cert = {"order": 2292, "construction": "spence", "group": [int(m) for m in S.G.mods],
        "sequences": [''.join('+' if x > 0 else '-' for x in s) for s in seqs], "P": Pm, "Q": H4}
Hm, _ = verify.build_spence(cert)
N = Hm.shape[0]
Hf = Hm.astype(np.float64)
G = Hf @ Hf.T
had = bool(np.isin(Hm, (-1, 1)).all() and np.array_equal(G, N * np.eye(N)))
print('H order', N, ' entries +-1 and H H^T = %d I:' % N, had)
if had:
    json.dump(cert, open('certificate_%s.json' % tag, 'w'))
    sha = hashlib.sha256()
    for ch in verify.csv_chunks(Hm): sha.update(ch)
    print('SHA-256 of the canonical CSV:', sha.hexdigest())

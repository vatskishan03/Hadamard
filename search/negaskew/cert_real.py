"""Certificate of a REAL-mode plug-in solution of qp (`real=1`: alphabet +-1, A real skew), a skew-Hadamard matrix of order 4v in the
cyclic Goethals-Seidel form of ../../skew/820 ('goethals-seidel', group Z_v, multiplier group H).
usage: python3 cert_real.py m q "[H_K generators]" "SOL line of qp real=1 (blocks A B C D, x0|x1|x2 with digits 0 and 2, one per K-orbit)"
The solution over Z_m x F_q is mapped to Z_v (v = m q) by the Chinese remainder theorem, the PAF sum and the skew relation are checked,
the matrix is built with the builder of ../../verify/verify.py (which also checks the complementarity), and out/cert_<4v>.json and
out/sequences_<4v>.txt are written (the certificate with its SHA-256).  A hit would be a new skew order 4v; none has been found.
Example (control, order 84):  ./qp shapes/z3f7.txt real=1 maxprint=1  (shape from  python3 qshape.py '[3]' 7 '[]' shapes/z3f7.txt)"""
import sys, json, os, hashlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', '..', 'verify'))
from qshape import kshape, cyc_gens
import verify as V
m = int(sys.argv[1]); q = int(sys.argv[2]); hK = json.loads(sys.argv[3]); sol = sys.argv[4]
sh = kshape([m], hK); v = m * q
blocks = sol.split(); blocks = blocks[1:] if blocks[0] == 'SOL' else blocks
squares = {x * x % q for x in range(1, q)}
seqs = []
for b in blocks:
    x0, x1, x2 = b.split('|'); x = []
    for g in range(v):
        o = sh['orb'][g % m]; f = g % q
        e = int(x0[o]) if f == 0 else (int(x1[o]) if f in squares else int(x2[o]))
        assert e in (0, 2); x.append(1 if e == 0 else -1)
    seqs.append(np.array(x))
gens = cyc_gens(m, q, hK)
H = {1}
while True:
    gr = H | {h * g % v for h in H for g in gens}
    if gr == H: break
    H = gr
H = sorted(H)
paf = [sum(int(np.dot(s, np.roll(s, -k))) for s in seqs) for k in range(v)]
A = seqs[0]
print('PAF peak', paf[0], 'max off-peak', max(abs(p) for p in paf[1:]), ' A skew', A[0] == 1 and all(A[(-g) % v] == -A[g] for g in range(1, v)))
reps = []
for s in seqs:
    seen = set(); r = []
    for g in range(v):
        if g in seen: continue
        orb = {g * h % v for h in H}; seen |= orb
        vals = {int(s[x]) for x in orb}; assert len(vals) == 1, 'not H-invariant'
        if s[g] == -1: r.append(g)
    reps.append(r)
cert = dict(order=4 * v, construction="goethals-seidel", skew=True, group=f"Z_{v}", v=v, multipliers=gens, subgroup_order=len(H), rowsums=[int(s.sum()) for s in seqs], blocks=reps)
Hm, _ = V.BUILDERS['goethals-seidel'](cert)      # raises ValueError unless the sequences are complementary
sha = hashlib.sha256()
for ch in V.csv_chunks(Hm): sha.update(ch)
cert['sha256'] = sha.hexdigest()
os.makedirs('out', exist_ok=True)
json.dump(cert, open(f'out/cert_{4 * v}.json', 'w'), indent=1)
open(f'out/sequences_{4 * v}.txt', 'w').write("\n".join("".join('+' if x > 0 else '-' for x in s) for s in seqs) + "\n")
n = 4 * v
print('order', n, 'sha256', cert['sha256'], '->', f'out/cert_{n}.json')
print('own check: +-1', bool(np.all(np.abs(Hm) == 1)), 'HH^T=nI', bool(np.array_equal(Hm @ Hm.T, n * np.eye(n, dtype=Hm.dtype))), 'skew', bool(np.array_equal(Hm + Hm.T, 2 * np.eye(n, dtype=Hm.dtype))))

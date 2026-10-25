"""Build the symmetric conference matrix of order 9p+1 from a bipcyc solution line and check it exactly.
usage: python3 build_bip.py p flip 'SOL g=.. R1: (d,cX,phY) ... R2: ...' [out.npy]"""
import re, sys
import numpy as np
p, flip, line = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
m = (p - 1) // 2
vec = [(a, b) for a in range(3) for b in range(3) if a or b]
neg = lambda v: ((-v[0]) % 3, (-v[1]) % 3)
posrep = {'A1': (0, 1), 'A2': (1, 0), 'B1': (1, 1), 'B2': (1, 2)}
def cls(v):
    for k, r in posrep.items():
        if v == r or v == neg(r): return k
R1s, R2s = line.split('R1:')[1].split('R2:')
blocks = {}
for rowname, txt in (('A1', R1s), ('A2', R2s)):
    a = posrep[rowname]
    for dd, col, ph in re.findall(r'\((\d+),c(\d),ph(\d)\)', txt):
        dd, col, ph = int(dd), int(col), int(ph)
        y, o = col >> 1, col & 1
        beta = posrep['B1' if y == 0 else 'B2']
        if o: beta = neg(beta)
        if dd <= m:                      # A at lower end: (u, v, c) = (a, -beta, -ph)
            s, u, v, c = dd, a, neg(beta), (-ph) % 3
        else:                            # A at upper end: (u, v, c) = (beta, -a, ph)
            s, u, v, c = p - dd, beta, neg(a), ph
        if u not in posrep.values():      # canonical: (u, v, c) ~ (-u, -v, -c)
            u, v, c = neg(u), neg(v), (-c) % 3
        assert s not in blocks
        blocks[s] = (u, v, c)
assert sorted(blocks) == list(range(1, m + 1)), sorted(blocks)
pts = [(a, b) for a in range(3) for b in range(3)]
def hyper(u, v, c):
    return np.array([[int((u[0]*x[0] + u[1]*x[1] + v[0]*y[0] + v[1]*y[1] - c) % 3 == 0) for y in pts] for x in pts], np.int64)
def T_dir(w):
    return np.array([[int(any(((k*w[0]) % 3, (k*w[1]) % 3) == ((x[0]-y[0]) % 3, (x[1]-y[1]) % 3) for k in range(3))) for y in pts] for x in pts], np.int64)
rook = T_dir((1, 1)) + T_dir((1, 2)) - 2 * np.eye(9, dtype=np.int64)
qr = {x * x % p for x in range(1, p)}
J9 = np.ones((9, 9), np.int64)
fsflip = 1 - flip                       # bipcyc sigma = +chi*flipsign; fs3 sigma = -chi*flipsign
sigf = lambda s: (-1 if s in qr else 1) * (-1 if fsflip else 1)
def Aof(s):
    u, v, c = blocks[s]; M = hyper(u, v, c)
    return M if sigf(s) == -1 else J9 - M
n = 9 * p + 1
C = np.zeros((n, n), np.int64)
for i in range(p):
    for j in range(p):
        dd = (j - i) % p
        if dd == 0: S = J9 - np.eye(9, dtype=np.int64) - 2 * rook
        elif dd <= m: S = J9 - 2 * Aof(dd)
        else: S = J9 - 2 * Aof(p - dd).T
        C[9*i:9*i+9, 9*j:9*j+9] = S
C[n-1, :n-1] = 1; C[:n-1, n-1] = 1
ok = (C == C.T).all() and (np.diag(C) == 0).all() and set(np.unique(C)) <= {-1, 0, 1} and np.array_equal(C @ C, (n - 1) * np.eye(n, dtype=np.int64))
print(f"p={p}: order {n}: symmetric conference matrix: {ok}")
if ok and len(sys.argv) > 4: np.save(sys.argv[4], C)

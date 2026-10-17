"""NG pair search: H-invariant antiperiodic +-1 sequences a, b of length n (H = <g> <= U(Z_2n)), meet in the middle on
the twisted power spectrum |a(zeta)|^2 at representatives of the H-orbits on odd exponents (zeta^n = -1).
    python3 ng_mitm.py n g [max_hits]"""
import sys, numpy as np
n, g = int(sys.argv[1]), int(sys.argv[2]); maxhits = int(sys.argv[3]) if len(sys.argv) > 3 else 3
m = 2 * n
def orbit(x, g, m):
    o = [x]; y = x * g % m
    while y != x: o.append(y); y = y * g % m
    return o
seen = set(); pairs = []
for x in range(m):
    if x in seen: continue
    o = orbit(x, g, m)
    if (x + n) % m in o: sys.exit(f"orbit of {x} closed under +n: H not admissible")
    o2 = orbit((x + n) % m, g, m)
    seen.update(o); seen.update(o2)
    pairs.append(o)
P = len(pairs)
# basis vectors: e_i = indicator of orbit O_i (+1) and O_i + n (-1), as a length-n sequence
B = np.zeros((P, n))
for i, o in enumerate(pairs):
    for x in o:
        if x < n: B[i, x] += 1
        else: B[i, x - n] -= 1
oseen = set(); reps = []
for k in range(1, m, 2):
    if k in oseen: continue
    oseen.update(orbit(k, g, m)); reps.append(k)
Z = np.exp(1j * np.pi * np.outer(np.arange(n), reps) / n)      # n x R : zeta_r^p
C = B @ Z                                                        # P x R spectra of basis vectors
R = len(reps)
print(f"n = {n}, |H| = {len(orbit(1, g, m))}, sign pairs {P}, spectral reps {R}", flush=True)
signs_all = None
def spectra(start, cnt):
    idx = np.arange(start, start + cnt)
    Sg = ((idx[:, None] >> np.arange(P)) & 1) * -2 + 1          # +-1, bit set -> -1
    return Sg, np.abs(Sg @ C) ** 2
N = 1 << P
key = lambda X: [tuple(r) for r in np.rint(X * 8).astype(np.int64)]
table = {}
CH = 1 << 18
for st in range(0, N, CH):
    Sg, Pw = spectra(st, min(CH, N - st))
    for i, k in enumerate(key(Pw)):
        table.setdefault(k, st + i)
hits = 0
for st in range(0, N, CH):
    Sg, Pw = spectra(st, min(CH, N - st))
    for i, k in enumerate(key(2 * n - Pw)):
        j = table.get(k)
        if j is None: continue
        sa = ((j >> np.arange(P)) & 1) * -2 + 1; sb = Sg[i]
        a = sa @ B; b = sb @ B
        # exact check in Z[x]/(x^n+1)
        def nacf(s):
            full = np.convolve(s, s[::-1]).astype(np.int64)   # index n-1 is shift 0
            r = full[n - 1:].copy()
            r[1:] -= full[:n - 1][::-1]
            return r
        c = nacf(a.astype(np.int64)) + nacf(b.astype(np.int64))
        if c[0] == 2 * n and not c[1:].any():
            hits += 1
            print("NG pair:", "".join('+' if v > 0 else '-' for v in a), "".join('+' if v > 0 else '-' for v in b), flush=True)
            if hits >= maxhits: sys.exit(0)
print("hits", hits, flush=True)

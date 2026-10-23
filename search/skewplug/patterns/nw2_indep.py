"""INDEPENDENT numpy check (no gl1.c) of the level-1 emptiness of one row-sum class of the N_w = 2 shapes K = Z_m, q = 2m + 1, H_K = 1, allocation (0,0,1,1), m-split (0,0,0,1):
the sum of the four per-block lag vectors F_A + F_B + F_C + F_D of PAF_u (u = s + q e) over the (m-1)/2 lag classes must vanish; hash of the vectors by a random linear form, sorted join A+B against C+D.
The explicit lists of the class are enumerated here from scratch (A skew, B without defect, C with a p-type defect at 0, D with an m-type defect at 0).
usage (from search/skewplug):  python3 patterns/nw2_indep.py PATTERN_INDEX [m q]       default m, q = 23, 47 (v = 1081); PATTERN_INDEX = index into row_patterns_skew(sh, (0,0,1,1), (0,0,0,1))
Control: python3 patterns/nw2_indep.py 2 11 23 and ... 3 11 23 print 35640 and 16720 key matches (the counts of nw_scan.py 11 23); for 1081 every class prints 'level-1 key matches 0'."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
import time, itertools, math, warnings
warnings.filterwarnings('ignore')
import numpy as np
from sk import SShape
from sk_pat import row_patterns_skew
pi = int(sys.argv[1])
N, q = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (23, 47)
sh = SShape([N], q, []); pats = row_patterns_skew(sh, (0, 0, 1, 1), (0, 0, 0, 1)); pat = pats[pi]
lags = [k for k in sh.repsPM if k != 0]; SHs = [np.array(sh.shift[k]) for k in lags]
def seq_list(n, total):
    k = (n + total) // 2
    out = np.full((math.comb(n, k), n), -1, dtype=np.int8)
    for r, comb in enumerate(itertools.combinations(range(n), k)): out[r, list(comb)] = 1
    return out
def feats(U):
    U = U.astype(np.int32); out = np.zeros((len(U), len(lags)), dtype=np.int32)
    for j, sh_ in enumerate(SHs): out[:, j] = (U * U[:, sh_]).sum(1)
    return out
def block(sig, eps, typ):
    if typ is None: return seq_list(N, sig).astype(np.int32)
    e0 = eps; s0 = 0 if typ == 0 else -2 * e0
    R = seq_list(N - 1, sig - s0); U = np.zeros((len(R), N), dtype=np.int32); U[:, 1:] = R; U[:, 0] = s0 + q * e0
    return U
(sB, eB), (sC, eC), (sD, eD) = pat
UB = block(sB, eB, None); UC = block(sC, eC, 0); UD = block(sD, eD, 1)
# A: skew (s on pairs), s(0) = 1
pairs = [(k, N - k) for k in range(1, (N - 1) // 2 + 1)]
ids = np.arange(1 << len(pairs)); UA = np.zeros((len(ids), N), dtype=np.int32); UA[:, 0] = 1
for j, (a, b) in enumerate(pairs):
    sg = 1 - 2 * ((ids >> j) & 1); UA[:, a] = sg; UA[:, b] = -sg
FA, FB, FC, FD = feats(UA), feats(UB), feats(UC), feats(UD)
print('sizes', len(FA), len(FB), len(FC), len(FD), flush=True)
rng = np.random.default_rng(3); mult = rng.integers(1, 2 ** 61, FA.shape[1], dtype=np.int64) | 1
key = lambda F: (F.astype(np.int64) * mult).sum(1)
# (A + B) table
kab = (key(FA)[:, None] + key(FB)[None, :]).ravel(); kab.sort()
kc, kd = key(FC), key(FD)
t0 = time.time(); hits = 0
step = max(1, 20_000_000 // len(kd))
for c0 in range(0, len(kc), step):
    tot = -(kc[c0:c0 + step, None] + kd[None, :]).ravel()
    hits += int((np.searchsorted(kab, tot, 'right') - np.searchsorted(kab, tot, 'left')).sum())      # all (a, b) matching each (c, d)
print(f'pattern {pi}: level-1 key matches {hits} ({time.time()-t0:.0f}s)', flush=True)

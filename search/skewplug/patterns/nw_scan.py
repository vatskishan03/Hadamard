"""EXACT level-1 count for the H_K = 1 shapes with at most one defect per free block (N_w = 2, 3: allocations (0,0,1,1), (0,1,1,1)): for every exact row-sum pattern the
explicit class-restricted lists (row-sum class of the pattern, the defect at position 0 by translation) and the exact hash join of ../gl1.c.
usage (from search/skewplug):  python3 patterns/nw_scan.py m q [maxlist]      (K = Z_m, q = 2m + 1 for N_w = 2)
Examples: nw_scan.py 11 23 -> 'level-1 instances' 35640 and 16720 for the two classes (the m = 11 patterns that survive the N_w = 2 lemma; SUMMARY total 52360, a few seconds);
nw_scan.py 15 31 -> 1,110,720 and 460,800 (total 1,571,520);  nw_scan.py 23 47 -> every class has 0 level-1 instances (1081 = 23 * 47; the lists are large: A 2048, B up to 1.35 M, a few hours)."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
import time, itertools, math, warnings
warnings.filterwarnings('ignore')
import numpy as np
from sk import SShape
from sk_pat import row_patterns_skew, allocations, m_splits
from gl1 import L1, BlockSpec
N, q = int(sys.argv[1]), int(sys.argv[2]); maxlist = float(sys.argv[3]) if len(sys.argv) > 3 else 4e6
sh = SShape([N], q, []); l1 = L1(sh)
def seq_list(n, total):
    if (n + total) % 2 or abs(total) > n: return np.zeros((0, n), dtype=np.int8)
    k = (n + total) // 2
    out = np.full((math.comb(n, k), n), -1, dtype=np.int8)
    for r, comb in enumerate(itertools.combinations(range(n), k)): out[r, list(comb)] = 1
    return out
def csize(n, total):
    if (n + total) % 2 or abs(total) > n: return 0
    return math.comb(n, (n + total) // 2)
def block_size(sig, eps, typ):
    if typ is None: return csize(N, sig)
    s0 = 0 if typ == 0 else -2 * eps
    return csize(N - 1, sig - s0)
def block(sig, eps, typ):
    if typ is None:
        Sm = seq_list(N, sig); return Sm, np.zeros_like(Sm)
    e0 = eps; s0 = 0 if typ == 0 else -2 * e0
    R = seq_list(N - 1, sig - s0)
    Sm = np.zeros((len(R), N), dtype=np.int8); Em = np.zeros((len(R), N), dtype=np.int8)
    Sm[:, 1:] = R; Sm[:, 0] = s0; Em[:, 0] = e0
    return Sm, Em
total = 0; undecided = 0; npat = 0
for a, em in allocations(sh):
    if a[0] != 0 or max(a[1:]) > 1: continue
    for ms in m_splits(a, em):
        pats = row_patterns_skew(sh, a, ms)
        for pat in pats:
            npat += 1
            types = []
            for X in (1, 2, 3):
                types.append(None if a[X] == 0 else (1 if ms[X] == 1 else 0))
            sizes = [2 ** ((N - 1) // 2)] + [block_size(pat[i][0], pat[i][1], types[i]) for i in range(3)]
            order = sorted(range(4), key=lambda i: sizes[i])
            pairs = max(float(sizes[order[0]]) * sizes[order[2]], float(sizes[order[1]]) * sizes[order[3]])
            if sizes[order[3]] > maxlist or pairs > 2e11:
                print('alloc', a, 'ms', ms, 'pat', pat, 'sizes', sizes, 'SKIPPED (too large)', flush=True); undecided += 1; continue
            lists = [sh.enum_A(0, 0)] + [block(pat[i][0], pat[i][1], types[i]) for i in range(3)]
            wts = [0] + [a[X] for X in (1, 2, 3)]
            specs = [BlockSpec(wts[i], S=lists[i][0], E=lists[i][1], cansign=False) for i in order]
            t0 = time.time()
            cnt, _, _, _ = l1.run(specs, maxout=1, per_class=3e4)
            total += cnt
            print('alloc', a, 'ms', ms, 'pat', pat, 'sizes', sizes, 'level-1 instances', cnt, f'[{time.time()-t0:.0f}s]', flush=True)
print('SUMMARY', N, q, 'patterns', npat, 'total level-1 instances (canonical signs, classes run)', total, 'skipped classes', undecided)

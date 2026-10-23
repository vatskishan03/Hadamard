"""INDEPENDENT numpy count (no gl1.c) of the level-1 instances of one defect pattern: F_A + F_B + F_C + F_D = 0 over the lag classes (F = per-block PAF of u = s + q e), explicit lists with all signs
and canonical block signs (as pat_exh.solve_pattern), sorted hash join A+B against C+D with all multiplicities.  Compare with the n1 of run_pat.py for the same pattern.
usage (from search/skewplug):  python3 patterns/pat_indep.py 'K' q 'hunits' 'A' 'B' 'C' 'D'      (patterns in the representation of run_pat.py: A = [[pair index, type], ...], B, C, D = [[orbit, type], ...])
e.g. the pattern class 2 of 1209 (Z39 x F31, <16>) has 544 level-1 instances (class list printed by run_pat.py)"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
import json, time, warnings
warnings.filterwarnings('ignore')
import numpy as np
from sk import SShape
from pat_exh import free_list, A_list
if len(sys.argv) != 8:
    sys.exit(__doc__)
K = json.loads(sys.argv[1]); q = int(sys.argv[2]); hu = json.loads(sys.argv[3])
A = tuple(tuple(x) for x in json.loads(sys.argv[4])); B, C, D = [tuple(tuple(x) for x in json.loads(a)) for a in sys.argv[5:8]]
sh = SShape(K, q, hu); N = sh.N
lags = [k for k in sh.repsPM if k != 0]; SHs = [np.array(sh.shift[k]) for k in lags]
def feats(S, E):
    U = (S.astype(np.int32) + q * E.astype(np.int32)); out = np.zeros((len(U), len(lags)), dtype=np.int32)
    for j, sh_ in enumerate(SHs): out[:, j] = (U * U[:, sh_]).sum(1)
    return out
FA = feats(*A_list(sh, A)); FB = feats(*free_list(sh, B)); FC = feats(*free_list(sh, C)); FD = feats(*free_list(sh, D))
print('sizes', len(FA), len(FB), len(FC), len(FD), flush=True)
rng = np.random.default_rng(11); mult = rng.integers(1, 2 ** 61, FA.shape[1], dtype=np.int64) | 1
key = lambda F: (F.astype(np.int64) * mult).sum(1)
kab = (key(FA)[:, None] + key(FB)[None, :]).ravel(); kab.sort()
kc, kd = key(FC), key(FD); hits = 0; step = max(1, 20_000_000 // len(kd)); t0 = time.time()
for c0 in range(0, len(kc), step):
    tot = -(kc[c0:c0 + step, None] + kd[None, :]).ravel()
    idx = np.searchsorted(kab, tot); idx[idx >= len(kab)] = len(kab) - 1
    cand = tot[kab[idx] == tot]
    for t in cand:
        lo = np.searchsorted(kab, t, 'left'); hi = np.searchsorted(kab, t, 'right'); hits += int(hi - lo)
print('independent level-1 instance count:', hits, f'({time.time()-t0:.0f}s)', flush=True)

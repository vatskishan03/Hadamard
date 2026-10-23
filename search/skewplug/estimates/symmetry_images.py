"""The symmetries of the K-level system: all 128 images of a level-1 instance have level-2 completions.  The K-level equations are sums of per-block autocorrelation-type vectors, so each block
may be reflected (k -> -k; the G-automorphism (k, f) -> (-k, -f)) and B, C, D may also be negated, separately (2 x 2^3 x 2^3 = 128 images of an instance of H_K-orbit-constant blocks; translations by
points of K fixed by H_K, tested in trans_test.py, add a factor |Fix|^3).  The raw level-1 counts of the samplers therefore overstate the number of distinct instances (prod_dedup.py removes the copies).
Here the instance is that of the repository's skew-Hadamard matrix skew/1892 (K = Z43, q = 11, H_K = <6>): its sequences are converted to K-level blocks (s, t, e), checked exactly at K level, the (s, e)
parts form the level-1 instance, and each of its 128 images is given to BOTH level-2 engines (skf and sl2), which must report at least one completion for every image.
usage (from search/skewplug):  python3 estimates/symmetry_images.py      (about 20 seconds; prints 'images with level-2 completions: skf 128 / 128, sl2 128 / 128')"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
import itertools, warnings
warnings.filterwarnings('ignore')
import numpy as np
from sk import SShape
from skf_io import skf_header, run_skf
from sl2io import header_text, run_sl2

order, m, q = 1892, 43, 11
v = order // 4
sh = SShape([m], q, [6])
txt = open(os.path.join(REPO, 'skew', str(order), 'sequences.txt')).read().split()
seqs = [np.array([1 if ch == '+' else -1 for ch in line]) for line in txt if set(line) <= set('+-') and len(line) == v][:4]
sq = {x * x % q for x in range(1, q)}; nsq = next(x for x in range(2, q) if x not in sq)
crt = {(g % m, g % q): g for g in range(v)}                      # Z_v = Z_m x F_q
blocks = []
for sqn in seqs:
    x0 = np.array([sqn[crt[(k, 0)]] for k in range(m)]); x1 = np.array([sqn[crt[(k, 1)]] for k in range(m)]); x2 = np.array([sqn[crt[(k, nsq)]] for k in range(m)])
    e = (x1 + x2) // 2; t = (x1 - x2) // 2; s = x0 - e
    blocks.append((s, t, e))
print('K-level exact check of skew/%d:' % order, sh.check_solution(blocks), flush=True)
neg = np.asarray(sh.negperm, dtype=np.int64)
inst = [(np.asarray(b[0], dtype=np.int8), np.asarray(b[2], dtype=np.int8)) for b in blocks]
images = []
for ra, rb, rc, rd, nb, nc, nd in itertools.product((0, 1), repeat=7):
    out = []
    for X, (S, E) in enumerate(inst):
        refl = (ra, rb, rc, rd)[X]; sg = -1 if (0, nb, nc, nd)[X] else 1
        S2 = S[neg] if refl else S; E2 = E[neg] if refl else E
        out.append((np.ascontiguousarray(sg * S2).astype(np.int8), np.ascontiguousarray(sg * E2).astype(np.int8)))
    images.append(out)
keys = {tuple(tuple(S.tolist()) + tuple(E.tolist()) for S, E in im) for im in images}
print('images', len(images), 'distinct', len(keys), flush=True)
res1, _ = run_skf(sh, images, maxsol=1, tlimit=120.0, hdr=skf_header(sh))
res2, _ = run_sl2(sh, images, maxsol=1, tlimit=120.0, maxfib=1500000, stop_first=False, hdr=header_text(sh))
n1 = sum(1 for i in range(len(images)) if res1[i]['compl'] > 0); n2 = sum(1 for i in range(len(images)) if res2[i]['compl'] > 0)
same = all(res1[i]['compl'] == res2[i]['compl'] for i in range(len(images)))
print(f'images with level-2 completions: skf {n1} / {len(images)}, sl2 {n2} / {len(images)}; completion counts equal in both engines: {same}')

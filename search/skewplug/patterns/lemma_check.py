"""Checks the N_w = 2 complementarity lemma (nw2_lemma.py) and the K-level structure on a REAL solution of the repository: the skew-Hadamard matrix skew/<order> over Z_v = Z_m x F_q (cyclic, H = squares of F_q*)
is converted to K-level blocks (s, t, e), checked exactly at K level, and the defects are listed (position, type p/m, sign e, row sums).  For skew/3612 (v = 903 = 21 * 43, N_w = 2) block 1 has an
m-type defect (e = +1, the rest of s equal to +2) and block 2 a p-type defect (e = +1, rest -2): |rest_C| = |rest_D| and eps_C rest_C = -eps_D rest_D, as the lemma says.
usage (from search/skewplug):  python3 patterns/lemma_check.py [order m q]      default: 3612 21 43"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
import warnings
warnings.filterwarnings('ignore')
import numpy as np
from sk import SShape
order, m, q = [int(x) for x in sys.argv[1:4]] if len(sys.argv) > 3 else (3612, 21, 43); v = order // 4
txt = open(os.path.join(REPO, 'skew', str(order), 'sequences.txt')).read().split()
seqs = [np.array([1 if ch == '+' else -1 for ch in line]) for line in txt if set(line) <= set('+-') and len(line) == v][:4]
sq = {x * x % q for x in range(1, q)}; nsq = next(x for x in range(2, q) if x not in sq)
sh = SShape([m], q, [])
def crt(k, f):
    for g in range(v):
        if g % m == k and g % q == f: return g
gidx = [[crt(k, f) for f in range(q)] for k in range(m)]
blocks = []
for sqn in seqs:
    x0 = np.array([sqn[gidx[k][0]] for k in range(m)]); x1 = np.array([sqn[gidx[k][1]] for k in range(m)]); x2 = np.array([sqn[gidx[k][nsq]] for k in range(m)])
    e = (x1 + x2) // 2; t = (x1 - x2) // 2; s = x0 - e
    blocks.append((s, t, e))
print('K-level exact check:', sh.check_solution(blocks))
for X, (s, t, e) in enumerate(blocks):
    d = [int(k) for k in np.flatnonzero(e)]
    print('block', X, 'defects', d, 'types', ['m' if s[k] != 0 else 'p' for k in d], 'e', [int(e[k]) for k in d], 'sigma', int(s.sum()), 'rest (sum of s off defects)', int(s.sum() - sum(s[k] for k in d)))

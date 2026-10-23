"""The N_w = 2 lemma (EXACT).  For q = 2m + 1 the defect identity allows two w-positions only as allocation (0, 0, 1, 1): one p-type defect in C and one m-type defect in D.
Translating both to 0, the level-1 equation at a lag k != 0 reads q eps_C a_k + (q - 2) eps_D b_k + SS(k) = 0 with a_k = s_C(k) + s_C(-k), b_k likewise and
|SS(k)| <= 4m - 4 < 2(q - 2); hence a_k = 0 <=> b_k = 0 and eps_C a_k = -eps_D b_k, so |sigma_C| = |sigma_D + 2 eps_D| for the row sums sigma of s.  This program lists the exact
row-sum patterns of the allocation (0,0,1,1) and counts those that satisfy the lemma.  H_K = 1 throughout: with a larger H_K there are only fewer patterns.
usage (from search/skewplug):  python3 patterns/nw2_lemma.py [m,q ...]       default pairs (m, q): 11,23 15,31 23,47 29,59 35,71   (v = m q = 253, 465, 1081, 1711, 2485)
Expected: m = 11: 2 of 7 patterns survive, m = 15: 2 of 8, m = 23: 0 of 4, m = 29: 0 of 6, m = 35: 3 of 17 (so 1081 = 23 * 47 and 1711 = 29 * 59 have no solution in this product)."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
import warnings
warnings.filterwarnings('ignore')
from sk import SShape
from sk_pat import row_patterns_skew, allocations, m_splits
for (N, q) in [(int(a.split(",")[0]), int(a.split(",")[1])) for a in sys.argv[1:]] or [(11, 23), (15, 31), (23, 47), (29, 59), (35, 71)]:
    sh = SShape([N], q, [])
    assert q == 2 * N + 1 and 2 * (q - 2) > 4 * N - 4
    for a, em in allocations(sh):
        if a != (0, 0, 1, 1): continue
        for ms in m_splits(a, em):
            pats = row_patterns_skew(sh, a, ms)
            ok = [p for p in pats if abs(p[1][0]) == abs(p[2][0] + 2 * p[2][1])]
            print(f'm={N} q={q} alloc {a} ms {ms}: {len(pats)} row-sum patterns, {len(ok)} satisfy the lemma |sigma_C| = |sigma_D + 2 eps_D|', ok[:3], flush=True)

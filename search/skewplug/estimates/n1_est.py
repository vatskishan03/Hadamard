"""Unbiased estimate of the number N1 of level-1 instances of a slot from sampled lists: (instances found) / (sampling fraction of the lists of B, C, D), all row-sum patterns, K-symmetric shapes
(exact configuration counts est_skew.T_free).  Compare with the Gaussian value of est_skew.estimate (logN1): at the hit slot of v = 473 the two agree within a factor 1.5 (measured 9.0e6, relative sd 0.29,
against the Gaussian 1.3e7 = 10^7.12).
usage (from search/skewplug):  python3 estimates/n1_est.py 'K' q 'hunits' 'al' 'ms' LB LC LD ROUNDS SEED
e.g.  python3 estimates/n1_est.py "[43]" 11 "[[6]]" "[0,4,7,10]" "[0,4,6,10]" 10000 60000 60000 4 1       (about 1 minute; the Gaussian estimate: python3 estimates/est_slots.py ... or est_skew.estimate)"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
DATA = os.path.join(HERE, '..', 'data'); RUNS = os.path.join(HERE, '..', 'runs')
import json, random, math
import numpy as np
from sk import SShape
from sk_samp import level1
from gl1 import L1
import est_skew as es
K = json.loads(sys.argv[1]); q = int(sys.argv[2]); hu = json.loads(sys.argv[3]); al = tuple(json.loads(sys.argv[4])); ms = tuple(json.loads(sys.argv[5]))
LB, LC, LD = [int(float(x)) for x in sys.argv[6:9]]; rounds = int(sys.argv[9]); seed = int(sys.argv[10])
sh = SShape(K, q, hu); l1 = L1(sh)
A = sh.enum_A(al[0], ms[0]); TA = len(A[0])
TF = [es.T_free(sh, n, m, False) for n, m in zip(al[1:], ms[1:])]
rng = random.Random(seed); ests = []
for r in range(rounds):
    cnt, S, E, st = level1(sh, l1, al, ms, None, LB, LC, LD, seed=rng.randrange(10 ** 6), A_list=A, maxout=1000, hseed=r + 1)
    Ls = st[3:7]
    frac = (Ls[1] / TF[0]) * (Ls[2] / TF[1]) * (Ls[3] / TF[2])
    ests.append(cnt / frac if frac > 0 else 0.0)
    print(f'round {r}: instances {cnt}, list sizes {Ls}, fraction {frac:.3e}, N1 est {ests[-1]:.3e}', flush=True)
print('N1 estimate (ordered level-1 instances of the slot, all row-sum patterns): mean %.3e, rel.sd %.2f; Gaussian est_skew logN1 for comparison: see est files' % (np.mean(ests), np.std(ests) / max(np.mean(ests), 1e-300)), 'log10 = %.2f' % (math.log10(max(np.mean(ests), 1e-300))))

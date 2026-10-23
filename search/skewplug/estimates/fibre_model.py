"""Workload model of the exact level 2 for H_K = 1 (measured by sampling, HEURISTIC).  For a shape K = Z_m, q and an allocation (nA,nB,nC,nD) with m-split, per block: the number n_t of free t signs,
the fraction p of random t signs that pass the exact single-block PSD filter |s^(chi) + i sqrt(q) t^(chi)|^2 <= 4|G| (all chi != 0), the resulting fibre size 2^n_t * p, and the best two-pair join cost
min over pairings of max(|X||Y|, |Z||W|) (log2; the t-gauge halves the largest fibre).  Any exact four-list method must enumerate two pair lists, so this is the number of pair evaluations
per level-1 instance after the single-block filters ('measured fibres', NOTES §88: 2^44 at Z25 to 2^77 at Z43).
usage (from search/skewplug):  python3 estimates/fibre_model.py 'K' q [NAL | 'alloc' 'ms']
  NAL   number of random (allocation, m-split) slots of the identity list (seed 1, default 4);   'alloc' 'ms'  one explicit slot, e.g.  "[25]" 23 "[0,0,2,3]" "[0,0,2,3]"
Rows of the table of NOTES §88 (log2 best two-pair join cost): Z25/q=23 (575) alloc (0,0,2,3): 44.2;  Z29/23: 52.0;  Z31/19: 54.0;  Z33/19, Z33/23: 57.8 .. 59.6;  Z35/7 (245), alloc (2,0,6,22): 45.3;
Z35/23: 62.8;  Z41/23: 72.5;  Z43/23 (989): 77.4."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
DATA = os.path.join(HERE, '..', 'data'); RUNS = os.path.join(HERE, '..', 'runs')
import json, math, random, warnings
warnings.filterwarnings('ignore')
import numpy as np
from sk import SShape
from sk_pat import allocations, m_splits
import est_skew as es
K = json.loads(sys.argv[1]); q = int(sys.argv[2]); explicit = len(sys.argv) > 4
nal = int(sys.argv[3]) if len(sys.argv) > 3 and not explicit else 4
sh = SShape(K, q); N = sh.N; T = sh.T2
rng = np.random.default_rng(1); prng = random.Random(1)
pool = [(a, ms) for (a, em) in allocations(sh) for ms in m_splits(a, em)]
n = N; z = np.arange(n)
W = np.exp(2j * np.pi * np.outer(np.arange(1, n), z) / n)      # (n-1, n) characters
slots = [(tuple(json.loads(sys.argv[3])), tuple(json.loads(sys.argv[4])))] if explicit else prng.sample(pool, min(nal, len(pool)))
for (a, ms) in slots:
    blocks = [es.sample_A(sh, a[0], ms[0], 60, rng)] + [es.sample_free(sh, nn, mm, 60, rng) for nn, mm in zip(a[1:], ms[1:])]
    lf = []; ps = []
    for b, (S, Tt, E) in enumerate(blocks):
        passes = []; nts = []
        for i in range(len(S)):
            nd = (E[i] == 0)
            if b == 0:   # A: t even on pairs, t(0) free
                free = [0] + [int(sh.spos[o][0]) for (o, po) in sh.apairs if E[i, int(sh.spos[o][0])] == 0]
                nt = len(free)
                tv = rng.integers(0, 2, (50, nt)) * 2 - 1
                tfull = np.zeros((50, n))
                tfull[:, 0] = tv[:, 0]
                for j, (o, po) in enumerate([(o, po) for (o, po) in sh.apairs if E[i, int(sh.spos[o][0])] == 0]):
                    for p in sh.spos[o]: tfull[:, p] = tv[:, j + 1]
                    for p in sh.spos[po]: tfull[:, p] = tv[:, j + 1]
            else:
                nt = int(nd.sum()); tfull = np.zeros((50, n)); tv = rng.integers(0, 2, (50, nt)) * 2 - 1; tfull[:, nd] = tv
            sf = S[i].astype(float)
            w = (W @ sf)[None, :] + 1j * math.sqrt(q) * (tfull @ W.T)
            passes.append(float((np.abs(w) ** 2 <= T + 0.01).all(1).mean())); nts.append(nt)
        p = float(np.mean(passes)); nt = float(np.mean(nts)); lf.append(nt + math.log2(max(p, 1e-9))); ps.append(p)
    lf[int(np.argmax(lf))] -= 1.0      # t-gauge
    o = sorted(lf); c = max(o[0] + o[3], o[1] + o[2])
    print(f'K={K} q={q} alloc {a} ms {ms}: n_t per block {[round(x) for x in [lf[i] - math.log2(max(ps[i], 1e-9)) for i in range(4)]]}, PSD pass fractions {[round(p, 3) for p in ps]}, log2 fibre sizes {[round(x, 1) for x in lf]}, log2 best 2-pair join cost {c:.1f}', flush=True)

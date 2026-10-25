#!/usr/bin/env python3
"""estpools.py TAB C T n1:pool:f0 n2:pool:f0 n3:pool:f0 n4:pool:f0   (all roles rotated) -> Gaussian E estimates (full covariance and independent)"""
import sys, subprocess, os, math, numpy as np
import est
tab, C, T = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
here = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(1)
L = []; info = []
for k, sp in enumerate(sys.argv[4:8]):
    n, path, f0 = sp.split(':'); n = int(float(n))
    masks = np.fromfile(path, dtype=np.uint64)
    if n < len(masks): masks = masks[rng.permutation(len(masks))[:n]]
    masks.tofile(f'w/_e{k}.m')
    subprocess.run([os.path.join(here, 'flat2'), 'keys', tab, f0, f'w/_e{k}.m', f'w/_e{k}.k'], check=True, capture_output=True)
    K = np.fromfile(f'w/_e{k}.k', dtype=np.int16).reshape(-1, C).astype(float)
    L.append((K, True)); info.append((len(masks), K.std(), K.mean()))
for i in info: print("role necklaces %d  key std %.3f mean %.3f" % i)
lE, lsz, lpdf = est.log10E(L, T)
print("full-covariance: log10 E_tuples = %.2f (log10 tuples %.2f, log10 pdf %.2f)  E_classes = %.3g" % (lE, lsz, lpdf, 10 ** lE / C))
sg2 = sum(s * s for _, s, _ in info)
lind = lsz - (C - 1) * math.log10(math.sqrt(2 * math.pi) * math.sqrt(sg2))
print("independent-class approx: log10 E_tuples = %.2f   (full - indep = %.2f)" % (lind, lE - lind))

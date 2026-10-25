#!/usr/bin/env python3
"""Gaussian-lattice estimate of the number of solutions among given lists.
E = prod(list sizes) * N(T; sum mu_i, sum Sigma_i)  on C-1 coordinates (one class is implied by the row sums).
lists: for each role a (keys array [n][C] int, rotated flag) ; entries = n * (C if rotated else 1)."""
import numpy as np, math

def rot_all(K):
    n, C = K.shape
    return np.concatenate([np.roll(K, -k, axis=1) for k in range(C)])  # key[(e+k)%C]

def stats(K, rotated, drop=-1):
    C = K.shape[1]
    if rotated:
        # analytic mean/cov of the rotation-symmetrised list: mean_e = overall mean ; cov = circulant average
        mu = np.full(C, K.mean())
        Kc = K - mu
        # covariance circulant: cov(e,e') = mean over k, rows of Kc[x,(e+k)] Kc[x,(e'+k)]
        full = np.zeros((C, C))
        # c_d = mean_{x,f} Kc[x,f]Kc[x,f+d]
        cd = np.array([np.mean(Kc * np.roll(Kc, -d, axis=1)) for d in range(C)])
        for e in range(C):
            for e2 in range(C):
                full[e, e2] = cd[(e2 - e) % C]
        return mu, full, K.shape[0] * C
    mu = K.mean(0); Kc = K - mu
    return mu, Kc.T @ Kc / len(K), K.shape[0]

def log10E(lists, T, drop=None):
    """lists: [(K, rotated)] ; T: target per class (scalar)"""
    C = lists[0][0].shape[1]
    mu = np.zeros(C); S = np.zeros((C, C)); lsz = 0.0
    for K, rot in lists:
        m, c, n = stats(K, rot)
        mu += m; S += c; lsz += math.log10(n)
    idx = [e for e in range(C) if e != (C - 1 if drop is None else drop)]
    Si = S[np.ix_(idx, idx)]
    d = np.full(len(idx), float(T)) - mu[idx]
    sign, logdet = np.linalg.slogdet(Si)
    q = d @ np.linalg.solve(Si, d)
    logpdf = -0.5 * (len(idx) * math.log(2 * math.pi) + logdet + q)
    return lsz + logpdf / math.log(10), lsz, logpdf / math.log(10)

if __name__ == "__main__":
    import sys
    # usage: est.py C T f1:rot f2:rot f3:rot f4:rot  (int16 key files)
    C = int(sys.argv[1]); T = int(sys.argv[2])
    L = []
    for a in sys.argv[3:7]:
        f, r = a.split(':')
        L.append((np.fromfile(f, dtype=np.int16).reshape(-1, C).astype(float), r == '1'))
    print("log10 E = %.2f (log10 tuples %.2f, log10 pdf %.2f)" % log10E(L, T))

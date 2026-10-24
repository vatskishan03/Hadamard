"""est_plain.py -- Gaussian estimate (HEURISTIC) of the number of ordered 4-tuples of sequences of a plain cyclic Spence shape (all row sums +1) with sum PAF = -4 off 0.
Uses the exact list L (row sum +1), the empirical mean/covariance of the class-value vectors and the empirical lattice spacing per class (gcd of differences).
E_ordered(rho=+1) ~ L^4 * Cell * N(T ; 4 mu, 4 Sigma)."""
import sys, math
import numpy as np
import os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import plain_spence3 as ps
def est(v, H, maxL=2_000_000):
    lab, reps, O, size, cls = ps.structure(v, H)
    Ms = ps.forms(v, lab, reps, O, cls)
    X = ps.enum_list(O, size)
    L = len(X)
    if L == 0: return None
    if L > maxL:
        X = X[np.random.default_rng(1).choice(L, maxL, replace=False)]
    Xi = X.astype(np.int64)
    V = np.stack([np.einsum('ni,ij,nj->n', Xi, M, Xi) for M in Ms], 1)
    C = V.shape[1]
    mu = V.mean(0); Sg = np.cov(V.T) if C > 1 else np.array([[V.var()]])
    sp = []
    for c in range(C):
        u = np.unique(V[:, c]); d = np.diff(u)
        g = 0
        for x in d[:2000]: g = math.gcd(g, int(x))
        sp.append(g if g > 0 else 1)
    T = np.full(C, -4.0)
    S4 = 4 * Sg
    sign, logdet = np.linalg.slogdet(S4)
    d = T - 4 * mu
    quad = float(d @ np.linalg.solve(S4, d))
    logN = -0.5 * C * math.log(2 * math.pi) - 0.5 * logdet - 0.5 * quad
    logcell = sum(math.log(g) for g in sp)
    logE = 4 * math.log(L) + logcell + logN
    return dict(v=v, O=O, h=len(H), C=C, L=L, spacing=sp, log10E_ordered=logE / math.log(10), log10E_unordered=(logE - math.log(24)) / math.log(10), quad=quad)
if __name__ == '__main__':
    v = int(sys.argv[1]); Omin = int(sys.argv[2]); Omax = int(sys.argv[3])
    subs = ps.all_subgroups(v)
    out = []
    for H in subs:
        lab, reps, O, size, cls = ps.structure(v, H)
        if not (Omin <= O <= Omax): continue
        # cheap row-sum feasibility first
        tot = int(sum(size)); off = tot
        dp = np.zeros(2 * tot + 1); dp[off] = 1
        for s in size:
            s = int(s); new = np.zeros_like(dp); new[s:] += dp[:-s]; new[:-s] += dp[s:]; dp = new
        if dp[off + 1] == 0: continue
        r = est(v, H)
        if r: out.append(r); print({k: (round(x, 2) if isinstance(x, float) else x) for k, x in r.items()}, flush=True)

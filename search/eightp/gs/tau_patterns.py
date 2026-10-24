"""tau_patterns.py -- all exact row-pattern classes of a real plug-in shape: 4-tuples of block states (n, m, sigma, eps, tau) satisfying the identities of tau_filter.
usage: tau_patterns.py Kmods r hg   -> prints the allocations (n_i.m_i) with their row-sum patterns (sigma_i, eps_i, tau_i) up to block order and negation"""
import sys, os, math, json, itertools, collections
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gsk import *


def block_states_full(sh):
    """set of (n, m, sigma, eps, tau) for one block (pruned by global bounds); as numpy rows"""
    K = sh.kCD; r = sh.r
    w = [int(x) for x in K.w]; forced = [bool(x) for x in K.forced]
    Tw = [int(sum(K.sig[p] for p in K.pos[o])) for o in range(K.nu)]
    D = 4 * sh.N + int(sh.phi[0])
    Smax = int(math.isqrt(int(sh.T2))); Tmax = int(math.isqrt(int(sh.T2) // r)); Rmax = int(math.isqrt(int(sh.Rsq)))
    order = sorted(range(K.nu), key=lambda o: -w[o])
    cap_s = [0] * (K.nu + 1); cap_e = [0] * (K.nu + 1); cap_t = [0] * (K.nu + 1)
    for i in range(K.nu - 1, -1, -1):
        o = order[i]; cap_s[i] = cap_s[i + 1] + 2 * w[o]; cap_e[i] = cap_e[i + 1] + w[o]; cap_t[i] = cap_t[i + 1] + abs(Tw[o])
    N = sh.N
    R_ = [N + 1, N + 1, 2 * N + 2, N + 2, N + 2]       # offsets for n, m, s, e, t
    Wd = [2 * R_[0] + 1, 2 * R_[1] + 1, 2 * R_[2] + 1, 2 * R_[3] + 1, 2 * R_[4] + 1]
    def enc(n, m, s, e, t):
        k = (n.astype(np.int64) + R_[0])
        k = k * Wd[1] + (m + R_[1]); k = k * Wd[2] + (s + R_[2]); k = k * Wd[3] + (e + R_[3]); k = k * Wd[4] + (t + R_[4])
        return k
    def dec(k):
        t = k % Wd[4] - R_[4]; k = k // Wd[4]; e = k % Wd[3] - R_[3]; k = k // Wd[3]; s = k % Wd[2] - R_[2]; k = k // Wd[2]
        m = k % Wd[1] - R_[1]; n = k // Wd[1] - R_[0]
        return n, m, s, e, t
    states = np.array([enc(np.array([0]), 0, 0, 0, 0)[0]], dtype=np.int64)
    for i, o in enumerate(order):
        ell = w[o]; T = Tw[o]
        n, m, s, e, t = dec(states)
        outs = []
        if not forced[o]:
            for sv in (1, -1):
                for tv in (1, -1): outs.append(enc(n, m, s + sv * ell, e, t + tv * T))
        for ev in (1, -1):
            outs.append(enc(n + ell, m, s, e + ev * ell, t))
            outs.append(enc(n + ell, m + ell, s - 2 * ev * ell, e + ev * ell, t))
        st = np.unique(np.concatenate(outs))
        n, m, s, e, t = dec(st)
        c = (r + 1) * n - 4 * m
        ok = (c <= D) & (np.abs(s) <= Smax + cap_s[i + 1]) & (np.abs(t) <= Tmax + cap_t[i + 1]) & (np.abs(e) <= (Rmax + Smax + cap_s[i + 1]) // r + cap_e[i + 1] + 1)
        states = st[ok]
    n, m, s, e, t = dec(states)
    rho = s + r * e; c = (r + 1) * n - 4 * m
    ok = (np.abs(rho) <= Rmax) & (np.abs(s) <= Smax) & (np.abs(t) <= Tmax) & (c <= D)
    return np.unique(np.stack([n[ok], m[ok], s[ok], rho[ok], t[ok]], 1), axis=0)


def patterns(sh, maxpairs=2e8):
    r = sh.r; D = 4 * sh.N + int(sh.phi[0]); T2 = int(sh.T2); Rsq = int(sh.Rsq); eps = sh.eps
    S = block_states_full(sh)
    # features of a block: c, q = s^2 + r t^2, rho^2, s*t   ; remember the states behind each feature
    n, m, s, rho, t = S.T
    c = (r + 1) * n - 4 * m
    F = np.stack([c, s * s + r * t * t, rho * rho, (s * t) if eps == 1 else np.zeros_like(s)], 1)
    Fu, inv = np.unique(F, axis=0, return_inverse=True)
    inv = inv.ravel()
    nf = len(Fu)
    if nf * nf > maxpairs: return None, S
    tgt = np.array([D, T2, Rsq, 0], dtype=np.int64)
    def key(V): return (((V[:, 0] << 44) | (V[:, 1] << 29)) | (V[:, 2] << 21)) + (V[:, 3] + (1 << 20))
    I1, I2 = np.triu_indices(nf)
    P = Fu[I1] + Fu[I2]
    ok = (P[:, 0] <= D) & (P[:, 1] <= T2) & (P[:, 2] <= Rsq)
    P = P[ok]; I1 = I1[ok]; I2 = I2[ok]
    kP = key(P)
    o = np.argsort(kP, kind='stable'); kPs = kP[o]
    need = key(tgt[None, :] - P) if len(P) else np.zeros(0, np.int64)
    lo = np.searchsorted(kPs, need, 'left'); hi = np.searchsorted(kPs, need, 'right')
    res = []
    for a in np.flatnonzero(hi > lo):
        for j in o[lo[a]:hi[a]]:
            if j < a: continue
            res.append((int(I1[a]), int(I2[a]), int(I1[j]), int(I2[j])))
    return (res, Fu, inv, S), S


if __name__ == '__main__':
    Kmods = json.loads(sys.argv[1]); r = int(sys.argv[2])
    hg = [(tuple(u), ph) for u, ph in json.loads(sys.argv[3])] if len(sys.argv) > 3 else []
    sh = GShape(Kmods, r, hg, 'spence', None)
    out, S = patterns(sh)
    if out is None: print('too large'); sys.exit()
    res, Fu, inv, S = out
    print('block states', len(S), 'feature classes', len(Fu), 'feature 4-tuples (sorted pairs)', len(res))
    # map back to allocations: for each feature class, the set of (n, m) allocations and (sigma, rho, tau)
    byf = collections.defaultdict(list)
    for idx, f in enumerate(inv): byf[int(f)].append(idx)
    allocs = collections.Counter()
    for (a, b, c_, d) in res:
        for fs in itertools.product(byf[a], byf[b], byf[c_], byf[d]):
            al = tuple(sorted((int(S[i][0]), int(S[i][1])) for i in fs))
            allocs[al] += 1
    print('distinct (n,m) allocations', len(allocs))
    for al, cnt in sorted(allocs.items())[:60]:
        print('  ', ' '.join(f'{n}.{m}' for n, m in al), 'patterns', cnt)

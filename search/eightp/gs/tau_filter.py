"""tau_filter.py -- EXACT necessary condition for the real GS-type Paley plug-ins (four free blocks, K-level system of gsk.py, search/NOTES.md section 86): the joint row-sum identities with the
EXACT t-sum tau, not only its parity/range.

Block X = s I + t Q + e J on K x F_r (kind 'CD' of pk.py).  Per H_K-orbit o (size ell = w[o], signed size T = sum of sig over the orbit):
   non-defect (orbit not forced): s = +-1, t = +-1 free, e = 0 :  c += 0,           sigma += s*ell,       eps += 0,        tau += t*T
   p-defect  (s = 0, e = +-1)                                   :  c += (r+1)*ell,   sigma += 0,           eps += e*ell,    tau += 0
   m-defect  (s = -2e, e = +-1)                                 :  c += (r-3)*ell,   sigma += -2e*ell,     eps += e*ell,    tau += 0
(forced orbits are defects only.)  The identities (proved in search/NOTES.md section 86, checked by brute force there and in t_identity.py):
   sum_X c_X = 4|K| + phi(0)                                  (defect identity J(0))
   sum_X (sigma_X^2 + r tau_X^2) = T2 = 4|G| - phi(0)         (I(0) summed: sum of the row sums of I)
   sum_X sigma_X tau_X = 0                                    (only if eps = +1, r = 1 mod 4: Q(0)-sum)
   sum_X (sigma_X + r eps_X)^2 = Rsq                          (sum of U)
A shape for which no 4-tuple of block states satisfies all of them has NO solution (EXACT).  The existing gsk.row_patterns restricts tau only to parity and range.
"""
import sys, math, itertools
import numpy as np
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gsk import *


def block_states(sh, prune=True):
    """all reachable final states (c, sigma, rho, tau) of one block (rho = sigma + r eps) as unique rows of an int64 array, pruned to the global bounds"""
    K = sh.kCD; r = sh.r
    w = [int(x) for x in K.w]; forced = [bool(x) for x in K.forced]
    Tw = [int(sum(K.sig[p] for p in K.pos[o])) if True else w[o] for o in range(K.nu)]
    D = 4 * sh.N + int(sh.phi[0])
    Smax = int(math.isqrt(int(sh.T2))); Tmax = int(math.isqrt(int(sh.T2) // r)); Rmax = int(math.isqrt(int(sh.Rsq)))
    order = sorted(range(K.nu), key=lambda o: -w[o])
    # remaining capacity for bounds
    cap_s = [0] * (K.nu + 1); cap_e = [0] * (K.nu + 1); cap_t = [0] * (K.nu + 1)
    for i in range(K.nu - 1, -1, -1):
        o = order[i]
        cap_s[i] = cap_s[i + 1] + 2 * w[o]; cap_e[i] = cap_e[i + 1] + w[o]; cap_t[i] = cap_t[i + 1] + abs(Tw[o])
    Eo = sh.N + 2; SO = 2 * sh.N + 2
    NS = 2 * SO + 1; NE = 2 * Eo + 1
    def enc(c, s, e, t): return (((c.astype(np.int64) * NS + (s + SO)) * NE + (e + Eo)) * NE + (t + Eo))
    def dec(k):
        t = k % NE - Eo; k = k // NE; e = k % NE - Eo; k = k // NE; s = k % NS - SO; c = k // NS
        return c, s, e, t
    states = np.array([enc(np.array([0]), 0, 0, 0)[0]], dtype=np.int64)
    for i, o in enumerate(order):
        ell = w[o]; T = Tw[o]
        c, s, e, t = dec(states)
        outs = []
        if not forced[o]:
            for sv in (1, -1):
                for tv in (1, -1):
                    outs.append(enc(c, s + sv * ell, e, t + tv * T))
        for ev in (1, -1):
            outs.append(enc(c + (r + 1) * ell, s, e + ev * ell, t))                       # p-defect
            outs.append(enc(c + (r - 3) * ell, s - 2 * ev * ell, e + ev * ell, t))        # m-defect
        st = np.unique(np.concatenate(outs))
        if prune:
            c, s, e, t = dec(st)
            ok = (c <= D) & (np.abs(s) <= Smax + cap_s[i + 1]) & (np.abs(t) <= Tmax + cap_t[i + 1]) & (np.abs(e) <= (Rmax + Smax + cap_s[i + 1]) // r + cap_e[i + 1] + 1)
            st = st[ok]
        states = st
    c, s, e, t = dec(states)
    rho = s + r * e
    ok = (np.abs(rho) <= Rmax) & (np.abs(s) <= Smax) & (np.abs(t) <= Tmax) & (c <= D)
    return np.unique(np.stack([c[ok], s[ok], rho[ok], t[ok]], 1), axis=0)


def tau_alive(sh, verbose=False, maxstates=4e7, maxbucket=6e7):
    """returns (alive: True/False/None (None = too large), info).  Blocks' feature vectors (c, q = s^2 + r t^2, rho^2, s t) are joined as pairs bucketed by the exact c-sum."""
    r = sh.r; D = 4 * sh.N + int(sh.phi[0]); T2 = int(sh.T2); Rsq = int(sh.Rsq); eps = sh.eps
    try:
        st = block_states(sh)
    except MemoryError:
        return None, 'DP memory'
    if verbose: print('block states', len(st), flush=True)
    if len(st) == 0: return False, 'no block state'
    c, s_, rho, t = st[:, 0], st[:, 1], st[:, 2], st[:, 3]
    F = np.stack([c, s_ * s_ + r * t * t, rho * rho, (s_ * t) if eps == 1 else np.zeros_like(s_)], 1)
    F = np.unique(F, axis=0)
    n = len(F)
    if verbose: print('distinct features', n, flush=True)
    cvals = np.unique(F[:, 0])
    groups = {int(cv): np.flatnonzero(F[:, 0] == cv) for cv in cvals}
    def key3(Q, Rr, Sg): return ((Q << 29) | (Rr << 21)) + (Sg + (1 << 20))
    def pair_keys(cs):
        """keys (q, rho^2, s t) of all unordered pairs with c-sum cs, as an array (empty if none)"""
        out = []
        for c1, g1 in groups.items():
            c2 = cs - c1
            if c2 < c1 or c2 not in groups: continue
            g2 = groups[c2]
            if c1 == c2: iu, ju = np.triu_indices(len(g1)); I = g1[iu]; J = g1[ju]
            else: I = np.repeat(g1, len(g2)); J = np.tile(g2, len(g1))
            Q = F[I, 1] + F[J, 1]; Rr = F[I, 2] + F[J, 2]; Sg = F[I, 3] + F[J, 3]
            ok = (Q <= T2) & (Rr <= Rsq)
            if ok.any(): out.append(key3(Q[ok], Rr[ok], Sg[ok]))
        return np.concatenate(out) if out else np.zeros(0, np.int64)
    def npairs(cs):
        tot = 0
        for c1, g1 in groups.items():
            c2 = cs - c1
            if c2 < c1 or c2 not in groups: continue
            tot += len(g1) * (len(g1) + 1) // 2 if c1 == c2 else len(g1) * len(groups[c2])
        return tot
    cs_all = sorted(set(a + b for a in groups for b in groups if a <= b))
    cs_set = set(cs_all)
    for cs in cs_all:
        cs2 = D - cs
        if cs2 < cs or cs2 not in cs_set: continue
        if max(npairs(cs), npairs(cs2)) > maxbucket: return None, f'bucket too large {npairs(cs)}x{npairs(cs2)}'
        A = pair_keys(cs)
        if len(A) == 0: continue
        B = A if cs2 == cs else pair_keys(cs2)
        if len(B) == 0: continue
        # complement of key (q, rho^2, st) is (T2 - q, Rsq - rho^2, -st)
        Aq = A >> 29; Ar = (A >> 21) & 255; As = (A & ((1 << 21) - 1)) - (1 << 20)
        need = key3(T2 - Aq, Rsq - Ar, -As)
        Bs = np.sort(B)
        pos = np.searchsorted(Bs, need); pos = np.clip(pos, 0, len(Bs) - 1)
        if (Bs[pos] == need).any(): return True, dict(blockstates=len(st), features=n)
    return False, dict(blockstates=len(st), features=n)


if __name__ == '__main__':
    import json
    Kmods = json.loads(sys.argv[1]); r = int(sys.argv[2])
    hg = [(tuple(u), ph) for u, ph in json.loads(sys.argv[3])] if len(sys.argv) > 3 else []
    sh = GShape(Kmods, r, hg, 'spence', None)
    print(tau_alive(sh, verbose=True))

#!/usr/bin/env python3
"""Semi-analytic cost model for flat-pool rotation-virtual k-trees (calibrated against measurements; see README.md).

Quantities (all per structure v, |U|, row-sum pattern):
  m = (v-1)/|U| cosets, C = m/2 classes, rotation group Z_C acts on keys by cyclic shift (rotation by C = reflection).
  pool_i(phi)   = comb(m, npos_i) * p(phi)           flat vectors with Phi <= phi  (p from `flat rem`, REM sampling)
  necklaces     = pool / m ;   entries per list L = C * necklaces   (all four roles rotated; each solution CLASS appears C times)
  sigma_kappa   = sqrt(v * kphi * phi) / 4          std of one key coordinate of a flat list member (kphi calibrated)
  E_tuples      = prod L_i * (2 pi)^(-(C-1)/2) (2 sigma_kappa)^-(C-1) * corr     expected solution tuples in the lists
  E_classes     = E_tuples / C
  join (4-list k-tree on s merge classes, memory-limited chunking q per list):
     pairs per central bucket P_b = L_c^2 * pp^s,  pp = 1/(sqrt(2 pi) sqrt(2) sigma_kappa)
     solutions per bucket      Sol_b = E_c * rho_s,  rho_s = shape * (1/(sqrt(2 pi) sigma_kappa))^s
     W = [c_pair * 2 P_b + c_scan * 4 L_c + c_idx * 4 L_c] / Sol_b          (cpu-seconds if c_* in seconds)
  generation: n necklaces cost n_max/g0 * ln(1/(1-n/n_max)) cpu-seconds per role (coupon-collector saturation).
"""
import math, re, json, sys
from math import comb, log10, log, sqrt, pi, exp

def load_rem(path):
    tab = []
    for line in open(path):
        mm = re.match(r"Phi<([0-9.]+)\s+(\d+)\s+p=([0-9.e+-]+)", line)
        if mm and float(mm.group(3)) > 0:
            tab.append((float(mm.group(1)), float(mm.group(3)), int(mm.group(2))))
    return tab

def logp(tab, phi):
    """log10 p(Phi<=phi) by log-linear interpolation on the REM table (hits >= 3 only), extrapolating with the fitted power law below."""
    pts = [(x, log10(p)) for x, p, c in tab if c >= 10]
    if phi <= pts[0][0]:
        # power law fitted on the first 6 points
        xs = [log(x) for x, _ in pts[:6]]; ys = [y for _, y in pts[:6]]
        n = len(xs); mx = sum(xs) / n; my = sum(ys) / n
        a = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
        return my + a * (log(phi) - mx)
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        if x0 <= phi <= x1:
            return y0 + (y1 - y0) * (phi - x0) / (x1 - x0)
    return pts[-1][1]

def chi2_cdf(x, k):
    # regularised lower incomplete gamma via series
    a = k / 2.0; x = x / 2.0
    if x <= 0: return 0.0
    term = 1.0 / a; s = term; n = 0
    while abs(term) > 1e-14 * abs(s) and n < 2000:
        n += 1; term *= x / (a + n); s += term
    return min(1.0, s * exp(-x + a * log(x) - math.lgamma(a)))

class Model:
    def __init__(self, v, ulen, npos, tab, kphi=0.95, corr_log10=0.0, g0=None, J=3e7, c_scan=2.0e-8, c_idx=1.5e-7,
                 mem_bytes=3.0e9, bytes_entry=20, bytes_pair=16, shape=0.5, smax=9):
        self.v, self.ulen = v, ulen
        self.m = (v - 1) // ulen; self.C = self.m // 2
        self.npos = npos; self.tab = tab; self.kphi = kphi; self.corr = corr_log10
        self.g0 = g0 or (lambda phi: 3000.0)       # distinct necklaces per cpu-second at low coverage
        self.cpair = 1.0 / J; self.c_scan = c_scan; self.c_idx = c_idx
        self.mem = mem_bytes; self.be = bytes_entry; self.bp = bytes_pair; self.shape = shape; self.smax = smax

    def necklaces_max(self, phi, i):
        return comb(self.m, self.npos[i]) * 10 ** logp(self.tab, phi) / self.m

    def sigma(self, phi):
        return sqrt(self.v * self.kphi * phi) / 4.0

    def logE_tuples(self, phi, n):           # n = necklaces per role (list of 4)
        C = self.C; sg = self.sigma(phi)
        lsum = sum(log10(C * x) for x in n)
        return lsum - (C - 1) * log10(sqrt(2 * pi) * 2 * sg) + self.corr

    def gen_time(self, phi, n):
        t = 0.0
        for i, x in enumerate(n):
            nm = self.necklaces_max(phi, i)
            f = min(x / nm, 0.999999)
            t += nm / self.g0(phi) * (-log(1 - f)) if f > 0 else 0.0
        return t

    def join_time(self, phi, n, qmax=1 << 20, mode='chunk'):
        """best (q, s) memory-feasible join time in cpu-seconds for per-role necklace counts n."""
        C = self.C; sg = self.sigma(phi); L = [C * x for x in n]
        E = 10 ** self.logE_tuples(phi, n)
        best = None
        pp = 1.0 / (sqrt(2 * pi) * sqrt(2) * sg); rr = 1.0 / (sqrt(2 * pi) * sg)
        q = 1
        while q <= qmax:
            Lc = [x / q for x in L]
            # 'chunk': arbitrary necklace chunks, q^4 sub-problems each holding E/q^4 solutions (penalty ~q^2 in time);
            # 'partition': sub-lists selected by the exact values of a few classes (every solution lies in one of ~q^3 sub-problems,
            # each holding E/q^3 solutions, penalty ~q) -- an optimistic bound, not implemented.
            Ec = E / q ** (4 if mode == 'chunk' else 3)
            mem_lists = sum(Lc) * self.be
            if mem_lists < 0.75 * self.mem:
                for s in range(3, self.smax + 1):
                    Pb = (Lc[0] * Lc[1]) * pp ** s + (Lc[2] * Lc[3]) * pp ** s
                    if Pb * self.bp > 0.3 * self.mem: continue
                    sol = Ec * self.shape * rr ** s
                    if sol <= 0: continue
                    Btot = 1.0 / min(sol, 1.0)                    # buckets needed for one expected solution
                    Nq = q ** (4 if mode == 'chunk' else 3)       # sub-problems
                    per_quad = max(1.0, Btot / Nq)                # buckets processed inside one sub-problem
                    cost_b = self.cpair * Pb + self.c_scan * sum(Lc) + self.c_idx * sum(Lc) / per_quad
                    t = cost_b * Btot
                    if best is None or t < best[0]: best = (t, q, s, Pb, sol)
            q *= 2
        return best

    def total(self, phi, f, mode='chunk'):
        """f: fraction of each flat pool taken; returns dict"""
        n = [f * self.necklaces_max(phi, i) for i in range(4)]
        logE = self.logE_tuples(phi, n)
        Ecl = 10 ** (logE - log10(self.C))
        g = self.gen_time(phi, n)
        j = self.join_time(phi, n, mode=mode)
        return dict(phi=phi, f=f, n=n, logE_tuples=logE, E_classes=Ecl, gen=g, join=j)

if __name__ == "__main__":
    pass

#!/usr/bin/env python3
"""predict.py : scan (Phi threshold, fraction of the flat pool used) for the cheapest flat-pool k-tree, from measured inputs.
usage: predict.py tg_<v>_<ulen>.json [--mem GB] [--corr C0] [--unlimited]
Returns best total cpu-seconds subject to E_classes(calibrated) >= EMIN (default 1.5)."""
import sys, json, math
from math import log10, log, comb
import model as M

J = 1.82e7; C_SCAN = 5.8e-9; C_IDX = 1.6e-7        # measured at v=487 (rj.c, 1 core, nice)
def corr_log10(ncons):                              # calibration of the Gaussian estimate (ratio found/predicted, partial-class matches, v=487,311,1103)
    return log10(0.223) + (ncons - 18) * log10(0.90) if ncons >= 18 else log10(0.223) * ncons / 18.0

def load(fn, mem_gb=3.0):
    d = json.load(open(fn)); m = d['m']; C = d['C']; v = d['v']; ulen = d['ulen']
    tab = [(x, p, c) for x, c, p in d['rem']]
    npos = [r['npos'] for r in d['roles']]
    gen = {float(k): val for k, val in d['gen'].items()}
    def g0(phi):
        xs = sorted(gen); 
        if phi <= xs[0]:
            # g scales roughly like p^0.9 below the measured range
            return gen[xs[0]] * 10 ** (0.9 * (M.logp(tab, phi) - M.logp(tab, xs[0])))
        for a, b in zip(xs, xs[1:]):
            if a <= phi <= b: return 10 ** (log10(gen[a]) + (log10(gen[b]) - log10(gen[a])) * (phi - a) / (b - a))
        return gen[xs[-1]]
    mod = M.Model(v, ulen, npos, tab, kphi=0.93, corr_log10=0.5 + corr_log10(C - 1), g0=g0, J=J, c_scan=C_SCAN, c_idx=C_IDX,
                  mem_bytes=mem_gb * 1.1e9)
    return d, mod

def scan(mod, emin=1.5, phis=None, fs=None, mode='chunk'):
    best = None
    phis = phis or [0.12 + 0.01 * i for i in range(0, 40)]
    fs = fs or ([10 ** (-3 + 0.25 * i) for i in range(0, 12)] + [0.5, 0.7, 0.9])
    for phi in phis:
        for f in fs:
            if f > 0.95: continue
            r = mod.total(phi, f, mode)
            if r['E_classes'] < emin or r['join'] is None: continue
            tot = r['gen'] + r['join'][0]
            if best is None or tot < best[0]: best = (tot, r)
    return best

if __name__ == "__main__":
    fn = sys.argv[1]; mem = 3.0
    d, mod = load(fn, mem)
    b = scan(mod)
    if b is None: print("no feasible point"); sys.exit()
    tot, r = b
    print(json.dumps(dict(v=d['v'], ulen=d['ulen'], C=d['C'], pattern=d['pattern'], phi=r['phi'], f=r['f'], necklaces=[int(x) for x in r['n']],
                          E_classes=r['E_classes'], gen_cpu_s=r['gen'], join_cpu_s=r['join'][0], q=r['join'][1], s=r['join'][2], total_cpu_s=tot)))

r"""Expected number of H-invariant c parts for a Theorem-7 template (search/odm_g.py), by a Gaussian model.

For each block, the c parts with the template's row sum form a list of size |L_b| (counted exactly from binomials over
the orbit-variable sizes); its (PAF, X) vectors at the pair-orbit representatives are sampled, and the sum of four
independent draws is modelled as a Gaussian per coordinate whose values live on the lattice of the observed residues
(gcd of differences).  The expected number of solutions is
    E = prod_b |L_b| * prod_k  g_k / (sigma_k sqrt(2 pi))      (coordinates treated as independent).
Calibration: the known hits (1716: t = 39, H = <16>, b on 3Z_39) should give E of order 1 or more.

Usage: python3 search/odm_estimate.py GROUP HGENS m TEMPLATE_SOURCE [--bsub d] [--hb GENS] [--b4w W] [--samples N]
       (the templates are taken from `odm_g.py ... --dry` output in the same configuration)
"""
import argparse
import itertools
import math
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gsgroup import Group, Orbits  # noqa: E402


def count_rowsum(sizes, target):
    """number of +-1 assignments of variables with these (signed) sizes and row sum target"""
    dp = {0: 1}
    for z in sizes:
        nd = {}
        for s, c in dp.items():
            nd[s + z] = nd.get(s + z, 0) + c
            nd[s - z] = nd.get(s - z, 0) + c
        dp = nd
    return dp.get(target, 0)


def sample_block(O, var, U, target, nsamp, rng):
    vsize, Dg, Cx, Lv = O.tables(var, 0, U)
    nv = len(var)
    vs = np.array(vsize)
    out = []
    tries = 0
    while len(out) < nsamp and tries < 200 * nsamp:
        tries += 1
        v = rng.choice([-1, 1], size=nv)
        # repair the row sum greedily
        rs = int(v @ vs)
        for _ in range(4 * nv):
            if rs == target:
                break
            j = rng.integers(nv)
            nrs = rs - 2 * v[j] * vs[j]
            if abs(nrs - target) < abs(rs - target):
                v[j] = -v[j]
                rs = nrs
        if rs != target:
            continue
        P = Dg.sum(axis=0) + np.einsum('i,j,jik->k', v, v, np.tril(np.ones((nv, nv)), -1)[:, :, None] * Cx)
        X = v @ Lv
        out.append(np.concatenate([P, X]))
    return np.array(out), count_rowsum([int(x) for x in vsize], target)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('group')
    ap.add_argument('hgens')
    ap.add_argument('m', type=int)
    ap.add_argument('--bsub', type=int, default=0)
    ap.add_argument('--hb', default='')
    ap.add_argument('--b4w', type=int, default=-1)
    ap.add_argument('--norm', default='')
    ap.add_argument('--samples', type=int, default=3000)
    ap.add_argument('--maxtpl', type=int, default=5)
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    cmd = [sys.executable, os.path.join(here, 'odm_g.py'), a.group, a.hgens, str(a.m), '--dry']
    for opt in ('bsub', 'hb', 'b4w', 'norm'):
        val = getattr(a, opt)
        if val not in ('', 0, -1):
            cmd += [f'--{opt}', str(val)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    tpls = []
    for ln in r.stdout.splitlines():
        if ln.startswith('TPL'):
            body = ln[4:]
            lst_end = body.index(']]') + 2
            tpl = eval(body[:lst_end])
            feas = eval(body[lst_end:].rsplit(' ', 1)[0])
            tpls.append((tpl, feas))
    print(r.stderr.strip().splitlines()[-3:], file=sys.stderr)
    G = Group(a.group)
    hg = [x for x in a.hgens.split('|') if x and x != '1'] or ['1' if len(G.comps) == 1 else None]
    H = G.closure([G.aut_perm(x) for x in hg if x])
    O = Orbits(G, H)
    rng = np.random.default_rng(1)
    t = G.t
    for ti, (tpl, feas) in enumerate(tpls[:a.maxtpl]):
        U = []
        for blk in tpl:
            u = np.zeros(t, dtype=np.int64)
            for g, v in blk:
                u[g] = v
            U.append(u)
        for rows in feas[:2]:
            samples, sizes = [], []
            for i in range(4):
                bt = 'skew' if i == 0 else 'plain'
                var = O.variables(bt, U[i])
                S, n = sample_block(O, var, U[i], 0 if i == 0 else rows[i - 1], a.samples, rng)
                samples.append(S)
                sizes.append(n)
            k = min(len(s) for s in samples)
            tot = sum(s[:k] for s in samples)
            var_ = np.var(tot, axis=0)
            logE = sum(math.log10(max(n, 1)) for n in sizes)
            for c in range(tot.shape[1]):
                vals = np.unique(np.concatenate([s[:, c] for s in samples]))
                diffs = np.diff(vals)
                g = int(np.gcd.reduce(diffs.astype(np.int64))) if len(diffs) else 1
                sig = math.sqrt(max(var_[c], 1e-9))
                logE += math.log10(max(g, 1) / (sig * math.sqrt(2 * math.pi))) if sig > 0 else 0
            print(f'template {ti} rows {rows}: lists {[f"{n:.2e}" for n in sizes]}  log10 E = {logE:.1f}')


if __name__ == '__main__':
    main()

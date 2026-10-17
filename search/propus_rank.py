"""Gaussian estimates for the propus structures of search/propus.py beyond the exhaustive range.

For v and a multiplier group H <= Z_v^*: A (and D, unless --freeD) constant on the orbits of <H, -1>, B on the orbits
of H; for every row pattern r_A^2 + 2 r_B^2 + r_D^2 = 4v the expected number of solutions of
PAF_A + 2 PAF_B + PAF_D = 0 off 0 is the product of the pools times the Gaussian density of the summed class values at
0 (full covariance of sampled members, lattice spacing from the sampled values, on the hyperplane fixed by the
row sums).

  python3 search/propus_rank.py v1,v2,.. [--mino 20] [--maxo 80] [--samples 20000]   -> $HADAMARD_WORK/propus_rank.jsonl
  python3 search/propus_rank.py m1,m2,.. --qborder q ...     (Z_q-bordered problem on Z_m, order 4(m + q))
"""
import argparse, json, math, sys
from math import gcd, log10
import numpy as np
from pathlib import Path as _P
sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
import propus as PR


def estimate(stA, stB, stD, rows, samples, rng, cls_sizes=None, target=None):
    """log10 of the expected number of solutions of PAF_A + 2 PAF_B + PAF_D = target (one value per class; 0 if None)."""
    keys, logs = [], 0.0
    for st, r, w in ((stA, rows[0], 1), (stB, rows[1], 2), (stD, rows[2], 1)):
        n = st.pool(r)
        if n <= 0: return -math.inf, 0
        logs += log10(n)
        keys.append(w * st.class_values(st.sample(r, samples, rng)))
    mu = sum(k.mean(0) for k in keys); Sig = sum(np.cov(k.T) for k in keys)
    w_, U = np.linalg.eigh(Sig); keep = w_ > 1e-6 * w_.max()
    d0 = (0 if target is None else np.asarray(target, np.float64)) - mu
    if (~keep).any() and np.abs(U[:, ~keep].T @ d0).max() > 1e-6 * max(1.0, float(np.abs(d0).max())):
        return -math.inf, int(keep.sum())         # the target is off the affine span of the class values
    z = U[:, keep].T @ d0
    ld = -0.5 * (keep.sum() * math.log(2 * math.pi) + np.log(w_[keep]).sum()) - 0.5 * float((z * z / w_[keep]).sum())
    # Lattice of the summed class values: spacing g_c per class (gcd of the sampled differences). The row sums fix
    # sum_c |c| val_c, so the covariance is singular in the direction w = (|c|); the density lives on the hyperplane
    # w^perp, where the lattice points have covolume prod(g_c) * |w| / gcd(g_c w_c). (Before October 14 the lattice
    # factor was dropped whenever the covariance was singular, i.e. always: the estimates were too low by about
    # C log10(4).)
    gs = []
    for c in range(Sig.shape[0]):
        g = 0
        for k in keys:
            u = np.unique(np.rint(k[:, c]).astype(np.int64))
            if len(u) > 1: g = gcd(g, int(np.gcd.reduce(np.diff(u))))
        gs.append(max(g, 1))
    gs = np.array(gs, np.int64)
    lat = float(np.log10(gs).sum())
    if (~keep).sum() == 1 and cls_sizes is not None:
        w = np.asarray(cls_sizes, np.float64)
        if abs(float(U[:, ~keep][:, 0] @ w) / np.linalg.norm(w)) > 0.999:
            gw = np.gcd.reduce((gs * np.asarray(cls_sizes, np.int64)).astype(np.int64))
            lat += log10(float(np.linalg.norm(w)) / float(gw))
    return logs + ld / math.log(10) + lat, int(keep.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('vs'); ap.add_argument('--mino', type=int, default=20); ap.add_argument('--maxo', type=int, default=80)
    ap.add_argument('--samples', type=int, default=20000); ap.add_argument('--out', default='propus_rank.jsonl')
    ap.add_argument('--qborder', type=int, default=0, help='rank the Z_q-bordered problem on Z_v (q | v; 1 = border of §27)')
    a = ap.parse_args()
    rng = np.random.default_rng(1)
    if a.qborder:
        return rank_bordered(a, rng)
    for v in [int(x) for x in a.vs.split(',')]:
        for S, gs in PR.subgroup_gens(v).items():
            Hm = set(S) | {(-h) % v for h in S}
            Ob = len({min(t * h % v for h in S) for t in range(v)})
            Oad = len({min(t * h % v for h in Hm) for t in range(v)})
            if max(Ob, Oad) < a.mino or max(Ob, Oad) > a.maxo: continue
            for symB, freeD in ((True, False), (False, False), (False, True)):
                if symB and (v - 1) in S: continue
                try:
                    st_ad, st_b = PR.structures(v, gs, symB)
                except AssertionError:
                    continue
                stA, stD = st_ad, (st_b if freeD else st_ad)
                for rows in PR.patterns(v):
                    if not freeD and rows[0] > rows[2]: continue
                    lE, rk = estimate(stA, st_b, stD, rows, a.samples, rng, st_ad.csize)
                    rec = dict(v=v, gens=list(gs), H=len(S), symB=symB, freeD=freeD, OA=int(stA.O), OB=int(st_b.O), OD=int(stD.O),
                               C=int(st_ad.C), rows=list(rows), logE=round(lE, 2), rank=rk,
                               pools=[round(log10(max(st.pool(r), 1)), 1) for st, r in ((stA, rows[0]), (st_b, rows[1]), (stD, rows[2]))])
                    print(json.dumps(rec), flush=True)
                    with open(W.WORK / a.out, 'a') as f: f.write(json.dumps(rec) + '\n')



def rank_bordered(a, rng):
    """Z_q-bordered propus on Z_m (m = each listed value): A on the orbits of <H, -1>, B and D on the orbits of H (H may
    contain -1); target -4q on the classes inside qZ_m and 0 off it; rows from propus.qborder_patterns(q)."""
    q = a.qborder
    for m in [int(x) for x in a.vs.split(',')]:
        if m % q: continue
        for S, gs in PR.subgroup_gens(m).items():
            Ob = len({min(t * h % m for h in S) for t in range(m)})
            if Ob < a.mino or Ob > a.maxo: continue
            st_ad, st_b = PR.structures(m, gs, False)
            tgt = np.array([-4 * q if g % q == 0 else 0 for g in st_ad.reps], np.float64)
            for rows in PR.qborder_patterns(q):
                if min(st_ad.pool(rows[0]), st_b.pool(rows[1]), st_b.pool(rows[2])) == 0: continue
                lE, rk = estimate(st_ad, st_b, st_b, rows, a.samples, rng, st_ad.csize, tgt)
                rec = dict(m=m, v=m + q, q=q, gens=list(gs), H=len(S), minus1=(m - 1) in S, OA=int(st_ad.O), OB=int(st_b.O),
                           C=int(st_ad.C), rows=list(rows), logE=round(lE, 2), rank=rk,
                           pools=[round(log10(max(st.pool(r), 1)), 1) for st, r in ((st_ad, rows[0]), (st_b, rows[1]), (st_b, rows[2]))])
                print(json.dumps(rec), flush=True)
                with open(W.WORK / a.out, 'a') as f: f.write(json.dumps(rec) + '\n')


if __name__ == '__main__':
    main()

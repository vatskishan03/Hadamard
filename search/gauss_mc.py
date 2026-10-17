#!/usr/bin/env python3
"""Hybrid estimate of the number of Szöllősi solutions in a multiplier structure.

The summed class vector of the four blocks is sampled from random combinations of uniform block samples. Classes are
ranked by how much their Monte Carlo hit rate at the target exceeds the Gaussian prediction (the compression classes
of composite N come first: their values are sums of a few large, repeated terms). They are taken into an exact Monte
Carlo estimate one by one while at least MINHITS samples hit all of them; the remaining classes get a Gaussian
density with the full covariance, conditioned on the Monte Carlo classes, times the covolume of their joint lattice.
For composite N the diagonal estimate of ssxx.gauss / bsz.gauss can be orders of magnitude low (search/NOTES.md).

    python3 search/gauss_mc.py N g [g2] --rows a,b,c,d [--bordered] [--Kb 20000] [--K 2000000]
"""
import argparse, math, sys
from math import log10, pi
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from bnw import Struct
from gauss_lat import hnf_basis

KINDS = ["L", "L", "A", "A"]


def estimate(st, rows, bordered, Kb=20000, K=2000000, seed=1, minhits=40, verbose=False, kinds=None):
    kinds = kinds or KINDS
    rng = np.random.default_rng(seed)
    tg = -4 if bordered else 0
    lN, blocks = 0.0, []
    for kind, r in zip(kinds, rows):
        n = st.counts(kind).get(r, 0)
        if bordered and r == 0:
            n //= 2
        if n == 0:
            return None
        lN += log10(n)
        blocks.append(np.rint(st.values(kind, st.sample(kind, r, Kb, rng))).astype(np.int32))
    C = st.C
    # summed class vectors from random combinations of the block samples
    Ssum = np.zeros((K, C), np.int32)
    for V in blocks:
        Ssum += V[rng.integers(0, Kb, K)]
    mu = Ssum.mean(0); Sg = np.cov(Ssum[:200000].T.astype(float))
    var = np.diag(Sg)
    t = np.full(C, tg)
    fixed = [c for c in range(C) if var[c] <= 1e-9]
    if any(abs(mu[c] - tg) > 1e-6 for c in fixed):
        return None
    free = [c for c in range(C) if var[c] > 1e-9]
    # spacing of each class (from differences)
    D = (Ssum[1:20001] - Ssum[0]).astype(np.int64)
    sp = {c: int(np.gcd.reduce(np.abs(D[:, c]))) or 1 for c in free}
    if any((tg - int(Ssum[0, c])) % sp[c] for c in free):
        return ("coset",)
    # lumpiness: variance / spacing^2 small relative to... rank classes by P(hit) marginal vs gaussian: use MC hit
    # frequency per class; lumpy = classes whose marginal hit frequency exceeds the Gaussian prediction the most
    hits_marg = (Ssum[:, free] == tg).mean(0)
    gpred = np.array([sp[c] / math.sqrt(2 * pi * var[c]) * math.exp(-(tg - mu[c]) ** 2 / (2 * var[c])) for c in free])
    ratio = (hits_marg + 1e-9) / gpred
    order = [free[i] for i in np.argsort(-ratio)]
    # greedily add classes to the MC set while the joint hit count stays >= minhits
    mask = np.ones(K, bool); mc = []
    for c in order:
        m2 = mask & (Ssum[:, c] == tg)
        if m2.sum() < minhits:
            break
        mask = m2; mc.append(c)
    pmc = mask.mean() if mc else 1.0
    rest = [c for c in free if c not in mc]
    # one class is implied by the row sums: drop it from the rest (if the MC set is not empty the dependency may
    # already be used; drop only when the rest has full rank issues)
    if rest:
        # check the linear dependency: weights = class sizes
        rest_drop = rest[1:] if len(mc) == 0 else rest
        Srr = Sg[np.ix_(rest_drop, rest_drop)]
        if mc:
            Srk = Sg[np.ix_(rest_drop, mc)]; Skk = Sg[np.ix_(mc, mc)]
            # conditional on the MC classes at the target
            sol = np.linalg.lstsq(Skk, (t[mc] - mu[mc]), rcond=None)[0]
            mcond = mu[rest_drop] + Srk @ sol
            Scond = Srr - Srk @ np.linalg.lstsq(Skk, Srk.T, rcond=None)[0]
        else:
            mcond = mu[rest_drop]; Scond = Srr
        w, U = np.linalg.eigh(Scond)
        ok = w > 1e-6 * w.max()
        if (~ok).sum() > 1:
            pass
        # restrict to the nonsingular subspace (the row-sum dependency removes one direction)
        Wn = U[:, ok]; wn = w[ok]
        d = Wn.T @ (t[rest_drop] - mcond)
        logdens = -0.5 * ok.sum() * math.log(2 * pi) - 0.5 * np.log(wn).sum() - 0.5 * float((d * d / wn).sum())
        # lattice covolume of the rest coordinates (restricted to the same dimension)
        Dr = [[int(x) for x in row] for row in D[:3000][:, rest_drop]]
        B, piv = hnf_basis(Dr, len(rest_drop))
        covol = 1
        for k, b in enumerate(B):
            covol *= b[piv[k]]
        covol = abs(covol)
        # if a direction was dropped as singular, the covolume of the full lattice overcounts by one spacing
        if (~ok).sum() == 1 and len(B) == len(rest_drop):
            # approximate: divide by the spacing along the dropped direction (use the smallest pivot)
            covol /= min(abs(B[k][piv[k]]) for k in range(len(B)))
        prest = logdens / math.log(10) + log10(covol)
    else:
        prest = 0.0
    tot = lN + log10(pmc) + prest
    if verbose:
        print(f"   MC classes {len(mc)} of {len(free)} (P = {pmc:.3g}), rest {len(rest)} (log10 P = {prest:.2f}), lN {lN:.2f}")
    return tot, len(mc), pmc, prest


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("N", type=int); ap.add_argument("gens", type=int, nargs="+")
    ap.add_argument("--rows", required=True); ap.add_argument("--bordered", action="store_true")
    ap.add_argument("--Kb", type=int, default=20000); ap.add_argument("--K", type=int, default=2000000)
    a = ap.parse_args()
    st = Struct(a.N, a.gens, plain=not a.bordered)
    rows = tuple(int(x) for x in a.rows.split(","))
    r = estimate(st, rows, a.bordered, a.Kb, a.K, verbose=True)
    print(f"N = {a.N} H = <{a.gens}> rows {rows}: {r if r is None or isinstance(r[0], str) else f'log10 E hybrid {r[0]:.2f}'}")


if __name__ == "__main__":
    main()

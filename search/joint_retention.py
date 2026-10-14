#!/usr/bin/env python3
"""Joint retention of the low-energy filters under the Gaussian model: (x_C, x_D) of a solution are jointly Gaussian
given x_A + x_B + x_C + x_D = t, and anti-correlated; lowe_sz multiplies the marginal retentions.

    python3 search/joint_retention.py N g --rows a,b,c,d [--bordered] --frac fC fD   (kept fractions of the C and D pools)
"""
import argparse, math, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from bnw import Struct
import lowe_sz

KINDS = ["L", "L", "A", "A"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("N", type=int); ap.add_argument("gens", type=int, nargs="+")
    ap.add_argument("--rows", required=True); ap.add_argument("--bordered", action="store_true")
    ap.add_argument("--frac", type=float, nargs=2, required=True); ap.add_argument("--K", type=int, default=20000)
    ap.add_argument("--M", type=int, default=400000)
    a = ap.parse_args()
    rows = tuple(int(x) for x in a.rows.split(","))
    st = Struct(a.N, a.gens, plain=not a.bordered)
    rng = np.random.default_rng(1); C = st.C; t = np.full(C, float(st.target))
    mus, sigs = [], []
    for k, r in zip(KINDS, rows):
        mu, sg, _ = lowe_sz.stats(st, k, r, a.K, rng); mus.append(mu); sigs.append(sg)
    Sig = sum(sigs); _, Sp = lowe_sz.psd_pinv_sqrt(Sig)
    taus, Smps, zs = {}, {}, {}
    for k, f in zip((2, 3), a.frac):
        Sm = Sig - sigs[k]; _, Smp = lowe_sz.psd_pinv_sqrt(Sm); zk = t - (sum(mus) - mus[k])
        V = st.values(KINDS[k], st.sample(KINDS[k], rows[k], a.K, rng))
        Ep = np.einsum("ij,jk,ik->i", V - zk, Smp, V - zk)
        taus[k] = float(np.quantile(Ep, f)); Smps[k] = Smp; zs[k] = zk
    # joint conditional law of (x_C, x_D) given the sum: x = mu + Sigma_blk S^+ (t - sum mu) + noise
    # sample all four blocks jointly: x_k = mu_k + L_k e_k, then project onto the constraint sum = t
    def sample(k, n):
        w, U = np.linalg.eigh(sigs[k]); w = np.clip(w, 0, None)
        return mus[k] + (rng.standard_normal((n, C)) * np.sqrt(w)) @ U.T
    X = [sample(k, a.M) for k in range(4)]
    tot = sum(X)
    resid = t - tot                           # distribute the residual with weights Sigma_k Sig^+
    Xc = [X[k] + (resid @ (Sp @ sigs[k]).T) for k in range(4)]     # x_k + Sigma_k Sig^+ (t - sum)
    inC = np.einsum("ij,jk,ik->i", Xc[2] - zs[2], Smps[2], Xc[2] - zs[2]) <= taus[2]
    inD = np.einsum("ij,jk,ik->i", Xc[3] - zs[3], Smps[3], Xc[3] - zs[3]) <= taus[3]
    pc, pd, pj = inC.mean(), inD.mean(), (inC & inD).mean()
    print(f"N = {a.N} <{a.gens}> rows {rows} kept fractions {a.frac}: retention C {pc:.4g}, D {pd:.4g}, product {pc * pd:.4g}, "
          f"joint {pj:.4g} (joint / product {pj / max(pc * pd, 1e-300):.3g})")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Why does the low-energy filter lose the solutions of a composite structure? Uses the complete lists and the exact
solution set cached by retention_test.py.

    python3 search/retention_diag.py N g --rows a,b,c,d
"""
import argparse, sys
from math import gcd
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import workdir as W
from bnw import Struct, enum_block
import lowe_sz

KINDS = ["L", "L", "A", "A"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("N", type=int); ap.add_argument("gens", type=int, nargs="+")
    ap.add_argument("--rows", required=True); ap.add_argument("--K", type=int, default=20000)
    a = ap.parse_args()
    rows = tuple(int(x) for x in a.rows.split(","))
    st = Struct(a.N, a.gens, plain=True)
    lists = []
    for k, (kind, r) in enumerate(zip(KINDS, rows)):
        V, M = enum_block(st, kind, r, 0, f"rt{k}")
        lists.append(V.astype(np.float64))
    cache = Path(f"{W.WORK}/lowe/sols_{a.N}_{'_'.join(map(str, a.gens))}_{'_'.join(map(str, rows))}.npy")
    S = np.load(cache)
    print(f"{len(S)} solutions; distinct members per block: {[len(set(S[:, k].tolist())) for k in range(4)]}")
    C = st.C
    order = [gcd(int(s), a.N) for s in st.reps]
    print("class reps", st.reps, "\nclass sizes", st.csize.tolist(), "\ngcd(rep, N)", order)
    t = np.zeros(C)
    rng = np.random.default_rng(1)
    mus, sigs = [], []
    for k, r in zip(KINDS, rows):
        mu, sg, _ = lowe_sz.stats(st, k, r, a.K, rng); mus.append(mu); sigs.append(sg)
    # exact means/covariances from the complete lists (the pools themselves)
    emus = [L.mean(0) for L in lists]; esigs = [np.cov(L.T) for L in lists]
    print("max |sampled mean - exact mean| / sd:", [round(float(np.max(np.abs(mus[k] - emus[k]) / np.sqrt(np.diag(esigs[k]) + 1e-9))), 3) for k in range(4)])
    Sig = sum(sigs)
    for k in (2, 3):
        Sm = Sig - sigs[k]; _, Smp = lowe_sz.psd_pinv_sqrt(Sm)
        zk = t - (sum(mus) - mus[k])
        X = lists[k]
        E = np.einsum("ij,jk,ik->i", X - zk, Smp, X - zk)
        Es = E[S[:, k]]
        q = np.searchsorted(np.sort(E), Es) / len(E)
        print(f"\nblock {'ABCD'[k]} (pool {len(X)}): solution members' energy quantiles: min {q.min():.3f} median {np.median(q):.3f} max {q.max():.3f}")
        sd = np.sqrt(np.diag(Sm))
        dev = (X[S[:, k]] - zk) / sd
        devall = (X[rng.integers(0, len(X), 20000)] - zk) / sd
        print("  per class: mean standardized deviation of solution members (random members):")
        for c in range(C):
            print(f"   class {c:2d} rep {st.reps[c]:3d} gcd {order[c]:3d}: sol mean {dev[:, c].mean():7.2f} sd {dev[:, c].std():5.2f} | all mean {devall[:, c].mean():6.2f} sd {devall[:, c].std():5.2f}")
        # alternative filters
        unit = [c for c in range(C) if order[c] == 1]
        for name, cl in (("unit classes only", unit), ("non-unit only", [c for c in range(C) if order[c] != 1])):
            if not cl: continue
            Ssub = Sm[np.ix_(cl, cl)]; _, P = lowe_sz.psd_pinv_sqrt(Ssub)
            E2 = np.einsum("ij,jk,ik->i", X[:, cl] - zk[cl], P, X[:, cl] - zk[cl])
            q2 = np.searchsorted(np.sort(E2), E2[S[:, k]]) / len(E2)
            print(f"  filter on {name}: solution quantiles min {q2.min():.3f} median {np.median(q2):.3f} max {q2.max():.3f}")


if __name__ == "__main__":
    main()

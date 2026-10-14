#!/usr/bin/env python3
"""Test of lowe_sz's retention model on a structure with a known exact solution set: enumerate all four blocks, join
exactly (all solutions), then keep the L lowest-energy members of each big block as lowe_sz would and count the
solutions that survive; compare with the model's retention.

    python3 search/retention_test.py N g [g2] --rows a,b,c,d [--bordered] --L 2e5 [--energy units] [--fibre]

--fibre keeps the L lowest-energy members among those whose compressions pass lowe_sz's fibre filters, with the
statistics of the allowed sub-pool (as lowe_sz.py does by default).
"""
import argparse, math, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import workdir as W
from bnw import Struct, enum_block, join4
import lowe_sz

KINDS = ["L", "L", "A", "A"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("N", type=int); ap.add_argument("gens", type=int, nargs="+")
    ap.add_argument("--rows", required=True); ap.add_argument("--bordered", action="store_true")
    ap.add_argument("--L", type=float, nargs="+", default=[2e5]); ap.add_argument("--K", type=int, default=20000)
    ap.add_argument("--energy", choices=["full", "units"], default="full")
    ap.add_argument("--fibre", action="store_true")
    a = ap.parse_args()
    rows = tuple(int(x) for x in a.rows.split(","))
    st = Struct(a.N, a.gens, plain=not a.bordered)
    lists, allowed = [], []
    filt = lowe_sz.fibre_filters(st, rows, a.bordered) if a.fibre else [[] for _ in range(4)]
    for k, (kind, r) in enumerate(zip(KINDS, rows)):
        V, M = enum_block(st, kind, r, 1 if (a.bordered and r == 0) else 0, f"rt{k}")
        lists.append(V.astype(np.int64))
        O = len(st.sizes(kind))
        Y = 1 - 2 * ((M[:, None] >> np.arange(O, dtype=np.uint64)) & np.uint64(1)).astype(np.int64)  # bit set: -1
        allowed.append(lowe_sz.fibre_ok(filt[k], Y))
    sizes = [len(x) for x in lists]
    order = min(([0, 1, 2, 3], [0, 2, 1, 3], [0, 3, 1, 2]),
                key=lambda o: sizes[o[0]] * sizes[o[1]] + sizes[o[2]] * sizes[o[3]])
    cache = Path(f"{W.WORK}/lowe/sols_{a.N}_{'_'.join(map(str, a.gens))}_{'_'.join(map(str, rows))}.npy")
    if cache.exists():
        S = np.load(cache); total = len(S); sols = S
    else:
        sols, total = join4([lists[k].astype(np.int16) for k in order], st.C, 10 ** 7, None, target=st.target)
        S = np.zeros((len(sols), 4), np.int64)
        for i, s in enumerate(sols):
            for pos, k in enumerate(order):
                S[i, k] = s[pos]
        cache.parent.mkdir(parents=True, exist_ok=True); np.save(cache, S)
    print(f"N = {a.N} H = <{a.gens}> rows {rows}: lists {sizes}, exact solutions {total} (rows read {len(sols)})")
    if a.fibre:
        print(f"  fibres at {sorted({f[0] for fk in filt for f in fk})}: allowed per block "
              f"{[round(float(x.mean()), 4) for x in allowed]}; solution members allowed "
              f"{[bool(allowed[k][S[:, k]].all()) if len(S) else None for k in range(4)]}")
    # lowe's model
    rng = np.random.default_rng(1); C = st.C; t = np.full(C, float(st.target))
    mus, sigs = [], []
    for k, (kd, r) in enumerate(zip(KINDS, rows)):
        Y, _ = lowe_sz.fibre_sample(st, kd, r, a.K, rng, filt[k])
        V = st.values(kd, Y); mus.append(V.mean(0)); sigs.append(np.cov(V.T))
    Sig = sum(sigs); _, Sig_pinv = lowe_sz.psd_pinv_sqrt(Sig)
    for L in a.L:
        keepsets, rets = [], []
        for k in range(4):
            if allowed[k].sum() <= L:
                keepsets.append(allowed[k].copy()); rets.append(1.0); continue
            Sm = Sig - sigs[k]
            U = lowe_sz.unit_classes(st) if a.energy == "units" else []
            if U:
                Smp = np.zeros((C, C)); Smp[np.ix_(U, U)] = lowe_sz.psd_pinv_sqrt(Sm[np.ix_(U, U)])[1]
            else:
                _, Smp = lowe_sz.psd_pinv_sqrt(Sm)
            zk = t - (sum(mus) - mus[k])
            X = lists[k].astype(np.float64)
            E = np.einsum("ij,jk,ik->i", X - zk, Smp, X - zk)
            E[~allowed[k]] = np.inf
            thr = np.partition(E, int(L) - 1)[int(L) - 1]
            keepsets.append(E <= thr)
            ck = mus[k] + sigs[k] @ Sig_pinv @ (t - sum(mus))
            Ck = sigs[k] - sigs[k] @ Sig_pinv @ sigs[k]
            w, U = np.linalg.eigh((Ck + Ck.T) / 2); w = np.clip(w, 0, None)
            Xs = ck + (rng.standard_normal((a.K, C)) * np.sqrt(w)) @ U.T
            Es = np.einsum("ij,jk,ik->i", Xs - zk, Smp, Xs - zk)
            rets.append(float(np.mean(Es <= thr)))
        surv = np.ones(len(S), bool)
        actual = []
        for k in range(4):
            inside = keepsets[k][S[:, k]]
            actual.append(float(inside.mean()) if len(S) else float("nan"))
            surv &= inside
        print(f"  L = {L:.3g}: model retention {[round(x, 4) for x in rets]}, actual (per block, marginal) "
              f"{[round(x, 4) for x in actual]}; solutions kept {int(surv.sum())} of {len(S)}, "
              f"model {len(S) * np.prod(rets):.1f}")


if __name__ == "__main__":
    main()

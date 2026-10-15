#!/usr/bin/env python3
"""Value multiplicity of a Szöllősi structure: members of a block with the same class values have the same PAF
conditions and the same energy, so the low-energy lists keep or drop them together and every solution comes with all
such combinations. Per block, the size-biased multiplicity n Σ_v p_v² (n the pool size, p_v the share of the
class-value vector v) is measured from collisions among M uniform members; the product over the four blocks predicts
the number of solutions per distinct value quadruple (within 10% on six structures with exact solution sets, search/
NOTES.md §41). The hybrid estimate divided by it estimates the distinct value quadruples.

    python3 search/value_mult.py N g[,g2] a,b,c,d [--bordered] [--M 1e6]
"""
import argparse, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bnw import Struct

KINDS = ["L", "L", "A", "A"]


def block_mult(st, kind, r, M, rng, bordered=False):
    """(n Σ p_v², n) for the block of kind and row sum r, or (None, n) when M samples show no collision"""
    n = st.counts(kind).get(r, 0)
    if bordered and r == 0:
        n //= 2
    Y = st.sample(kind, r, M, rng)
    V = np.rint(st.values(kind, Y)).astype(np.int64)
    _, cnt = np.unique(V, axis=0, return_counts=True)
    pairs = float((cnt * (cnt - 1) / 2).sum())
    if pairs == 0:
        return None, n
    return n * pairs / (M * (M - 1) / 2), n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("N", type=int); ap.add_argument("gens"); ap.add_argument("rows")
    ap.add_argument("--bordered", action="store_true"); ap.add_argument("--M", type=float, default=1e6)
    ap.add_argument("--seed", type=int, default=11)
    a = ap.parse_args()
    st = Struct(a.N, [int(x) for x in a.gens.split(",")], plain=not a.bordered)
    rows = [int(x) for x in a.rows.split(",")]
    rng = np.random.default_rng(a.seed)
    ms = []
    for k, (kd, r) in enumerate(zip(KINDS, rows)):
        if k == 3 and rows[3] == rows[2]:
            ms.append(ms[2])
            print("block D: the pool of C"); continue
        m, n = block_mult(st, kd, r, int(a.M), rng, a.bordered)
        ms.append(m)
        print(f"block {'ABCD'[k]} ({kd}, row sum {r}): pool {n:.3g}, multiplicity "
              f"{'no collision in the sample' if m is None else f'{m:.2f}'}", flush=True)
    if all(m is not None for m in ms):
        print(f"product {np.prod(ms):.1f}: solutions per distinct value quadruple")


if __name__ == "__main__":
    main()

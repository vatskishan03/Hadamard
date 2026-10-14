#!/usr/bin/env python3
"""Exhaustive multiplier-structure search for binary Legendre pairs of a given odd length l.

A Legendre pair is (a, b), +-1 sequences of length l with PAF_a(k) + PAF_b(k) = -2 for k != 0 (row sums +-1). With
both sequences constant on the orbits of a multiplier group H <= Z_l^*, the condition holds on the <H, -1>-classes,
and for each H the two complete lists (row sum 1; -a has the same PAF) are joined exactly (sort join). Paley-product
groups (H trivial mod one factor, the squares mod a prime q = 3 (mod 4)) are included with every other subgroup.

    python3 search/lp_mult.py l [l ...] [--maxo 34] [--log FILE]
"""
import argparse, json, os, struct, subprocess, sys, time
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import workdir as W
from bnw import Struct, subgroups, mask_to_y, paf


def enum(st, r):
    d = W.WORK / "lpm"; d.mkdir(exist_ok=True)
    base = d / f"l{st.N}_{'_'.join(map(str, st.gens))}_{r}"
    vals, masks = base.with_suffix(".v"), base.with_suffix(".m")
    if not vals.exists():
        inp = base.with_suffix(".in")
        sz = st.sizes("A"); S = st.S["A"]
        with open(inp, "wb") as f:
            f.write(struct.pack("4i", len(sz), st.C, r, 0))
            f.write(np.asarray(sz, np.int32).tobytes()); f.write(np.asarray(S, np.int32).tobytes())
        subprocess.run([W.bin("bnw_enum"), str(inp), str(vals) + ".tmp", str(masks) + ".tmp"], check=True,
                       capture_output=True)
        os.replace(str(vals) + ".tmp", vals); os.replace(str(masks) + ".tmp", masks)
    return np.fromfile(vals, np.int16).reshape(-1, st.C), np.fromfile(masks, np.uint64)


def search(l, maxo=34, log=None, maxlist=2e8):
    found = []
    subs = sorted(subgroups(l).items(), key=lambda kv: -len(kv[0]))
    for H, gens in subs:
        st = Struct(l, gens, plain=True)
        O = len(st.sizesA)
        if O > maxo:
            continue
        n = st.counts("A").get(1, 0)
        if n == 0 or n > maxlist:
            continue
        t0 = time.time()
        V, M = enum(st, 1)
        key = lambda X: np.ascontiguousarray(X.astype(np.int16)).view(np.dtype((np.void, 2 * st.C))).ravel()
        ka = key(V); order = np.argsort(ka); ka = ka[order]
        kb = key(-2 - V.astype(np.int32))
        pos = np.searchsorted(ka, kb); pos[pos >= len(ka)] = 0
        hit = np.flatnonzero(ka[pos] == kb)
        print(f"l = {l}, H = <{','.join(map(str, gens))}> (order {len(H)}), {O} orbits, {st.C} classes: list {len(V)}, "
              f"{len(hit)} matches ({time.time() - t0:.1f}s)", flush=True)
        for h in hit[:3]:
            a = st.expand("A", mask_to_y(M[order[pos[h]]], O)); b = st.expand("A", mask_to_y(M[h], O))
            P = paf(a) + paf(b)
            if P[0] == 2 * l and np.all(P[1:] == -2):
                rec = {"l": l, "gens": list(gens), "a": "".join("+" if v > 0 else "-" for v in a),
                       "b": "".join("+" if v > 0 else "-" for v in b)}
                found.append(rec)
                print(f"  Legendre pair of length {l} verified", flush=True)
                if log:
                    with open(log, "a") as f:
                        f.write(json.dumps(rec) + "\n")
        if found:
            return found
    return found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("l", type=int, nargs="+"); ap.add_argument("--maxo", type=int, default=34)
    ap.add_argument("--log")
    a = ap.parse_args()
    for l in a.l:
        f = search(l, a.maxo, a.log)
        print(f"l = {l}: {'FOUND' if f else 'nothing'}", flush=True)


if __name__ == "__main__":
    main()

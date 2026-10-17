#!/usr/bin/env python3
"""Re-rank lowe-feasible structures of the frontier scan with the hybrid estimate (gauss_mc) and the in-list
estimate of lowe_sz (with its retention model, scaled by the hybrid/diagonal ratio). Patterns that fail the exact
compression check (compress_check.py) are written with "ruled_out" and do not count towards --top.

    python3 search/hybrid_rank.py FRONTIER.jsonl OUT.jsonl [--minlogE -6] [--maxpool 6e10] [--top 3]
"""
import argparse, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from bnw import Struct
import gauss_mc
import compress_check


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inp"); ap.add_argument("out")
    ap.add_argument("--minlogE", type=float, default=-6); ap.add_argument("--maxpool", type=float, default=6e10)
    ap.add_argument("--top", type=int, default=3); ap.add_argument("--skip", default="")
    a = ap.parse_args()
    skip = {int(x) for x in a.skip.split(",") if x}
    rows = [json.loads(l) for l in open(a.inp)]
    sel = [r for r in rows if r["work"] > 2e10 and max(r["sizes"][2:]) <= a.maxpool and r["logE"] >= a.minlogE
           and r["N"] not in skip]
    groups = {}
    for r in sel:
        groups.setdefault((r["N"], r["bordered"], tuple(r["gens"])), []).append(r)
    done = set()
    if Path(a.out).exists():
        for l in open(a.out):
            d = json.loads(l); done.add((d["N"], d["bordered"], tuple(d["gens"]), tuple(d["rows"])))
    order = sorted(groups.items(), key=lambda kv: -max(r["logE"] for r in kv[1]))
    with open(a.out, "a") as f:
        for (N, b, gens), rs in order:
            st = Struct(N, gens, plain=not b)
            n = 0
            for r in sorted(rs, key=lambda r: -r["logE"]):
                if n >= a.top:
                    break
                key = (N, b, tuple(gens), tuple(r["rows"]))
                if key in done:
                    n += 1; continue
                cc = compress_check.check("bsz" if b else "ssxx", N, list(gens), tuple(r["rows"]), 3_000_000, 64,
                                          stop=True)
                dead = [dd for dd, c in cc.items() if c == 0]
                if dead:
                    d = dict(r); d["hybrid"] = None; d["ruled_out"] = dead[0]
                    f.write(json.dumps(d) + "\n"); f.flush()
                    continue
                n += 1
                h = gauss_mc.estimate(st, tuple(r["rows"]), b, Kb=10000, K=1000000)
                hv = None if (h is None or isinstance(h[0], str)) else round(float(h[0]), 2)
                d = dict(r); d["hybrid"] = hv
                f.write(json.dumps(d) + "\n"); f.flush()
                print(N, "b" if b else "c", gens, r["rows"], "diag", r["logE"], "hybrid", hv,
                      ["%.2g" % s for s in r["sizes"]], flush=True)


if __name__ == "__main__":
    main()

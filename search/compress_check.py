#!/usr/bin/env python3
"""Necessary conditions from the compressions: for every prime-power divisor d of the group order, join the four
blocks' compressions to Z_d exactly (search/cchain4.py, level 1). An empty join rules the structure/pattern out.
A level whose state space exceeds maxstates is skipped (None); --dmax skips the levels above it (the large prime
levels take minutes and gigabytes, the small ones milliseconds).

    python3 search/compress_check.py ssxx N g[,g2] a,b,c,d [--dmax D] [--joint]  (cyclic Szöllősi)
    python3 search/compress_check.py bsz N g[,g2] r1,r2,r3,r4 [--dmax D]  (bordered, rows like 0,0,2,0)
    python3 search/compress_check.py --jsonl IN --out OUT [--dmax D]       (records with N, gens, rows, bordered,
                                                                            e.g. from frontier_scan.py; adds
                                                                            "compress" and "ruled_out")
"""
import argparse, json, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import cchain4 as CC
import chain as CH


def joint_level(mode, N, gens, pat, ds, maxstates=2_000_000, maxpairs=2e7):
    """joint check of several coprime levels: each block's joint compression states (to every Z_d of ds at once) are
    joined on the PAF conditions of all levels together. Stronger than the levels one by one: a member must realise
    all its compressions at once (at N = 51, H = <38>, it rules out three patterns that pass d = 3 and d = 17 alone).
    Number of joint solutions, 0, None (states too many) or "large"."""
    kinds, T, bordered = CC.MODES[mode]
    v = N - 1 if bordered else N
    ch = CC.KChain(v, list(gens), list(ds), T, kinds)
    S = [ch.states(k, r, maxstates) for k, r in zip(kinds, pat)]
    if any(x is None for x in S):
        return None
    U = [np.unique(K_, axis=0).astype(np.int64) for K_, _ in S]
    if any(len(u) == 0 for u in U):
        return 0
    def paf_all(u):
        return np.concatenate([ch.paf(lev, u[:, ch.splits[lev]:ch.splits[lev + 1]]) for lev in range(len(ch.ds))], axis=1)
    tgt = np.concatenate([ch.target(lev) for lev in range(len(ch.ds))])
    try:
        return len(CH.join4([paf_all(u) for u in U], tgt, maxpairs=maxpairs))
    except CH.JoinTooLarge:
        return "large"


def check(mode, N, gens, pat, maxstates=2_000_000, dmax=None, stop=False, joint=False, maxpairs=2e7):
    """{d: number of level-1 solutions, 0 (ruled out), None (state space too large) or "large" (join too large)};
    stop: return at the first empty level. joint: when every level passes, also join each pair of coprime levels
    that could be listed (key "d1x d2")"""
    kinds, T, bordered = CC.MODES[mode]
    v = N - 1 if bordered else N
    ps = [p for p in range(2, v) if v % p == 0 and all(p % r for r in range(2, int(p ** 0.5) + 1))]
    cand = sorted({p ** e for p in ps for e in range(1, 8) if v % (p ** e) == 0 and p ** e < v})
    out = {}
    for d in cand:
        if dmax and d > dmax:
            continue
        c1 = CC.KChain(v, gens, [d], T, kinds)
        Sd = [c1.states(k, r, maxstates) for k, r in zip(kinds, pat)]
        if any(x is None for x in Sd):
            out[d] = None; continue
        U = [np.unique(K_, axis=0).astype(np.int64) for K_, _ in Sd]
        try:
            out[d] = 0 if any(len(u) == 0 for u in U) else len(CH.join4([c1.paf(0, u) for u in U], c1.target(0)))
        except CH.JoinTooLarge:
            out[d] = "large"
        if stop and out[d] == 0:
            break
    if joint and not any(c == 0 for c in out.values()):
        ok = [d for d, c in out.items() if isinstance(c, int) and c > 0]
        for i in range(len(ok)):
            for j in range(i + 1, len(ok)):
                if np.gcd(ok[i], ok[j]) != 1:
                    continue
                key = f"{ok[i]}x{ok[j]}"
                out[key] = joint_level(mode, N, gens, pat, [ok[i], ok[j]], maxstates, maxpairs)
                if stop and out[key] == 0:
                    return out
    return out


def screen(path, out, maxstates, dmax, joint=False):
    n = bad = 0
    with open(out, "a") as fo:
        for line in open(path):
            r = json.loads(line)
            pat = tuple(r.get("rows") or r.get("pat"))
            t0 = time.time()
            res = check("bsz" if r.get("bordered") else "ssxx", r["N"], list(r["gens"]), pat, maxstates, dmax, stop=True,
                        joint=joint)
            r["compress"] = {str(d): c for d, c in res.items()}
            r["ruled_out"] = next((d for d, c in res.items() if c == 0), None)
            r["compress_seconds"] = round(time.time() - t0, 1)
            fo.write(json.dumps(r) + "\n"); fo.flush()
            n += 1; bad += r["ruled_out"] is not None
    print(f"{n} records, {bad} ruled out")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", nargs="?", choices=["ssxx", "bsz"]); ap.add_argument("N", nargs="?", type=int)
    ap.add_argument("gens", nargs="?"); ap.add_argument("pat", nargs="?")
    ap.add_argument("--dmax", type=int); ap.add_argument("--maxstates", type=float, default=2e6)
    ap.add_argument("--jsonl"); ap.add_argument("--out")
    ap.add_argument("--joint", action="store_true", help="also join pairs of coprime levels jointly")
    a = ap.parse_args()
    if a.jsonl:
        screen(a.jsonl, a.out, int(a.maxstates), a.dmax, a.joint); return
    mode, N, gens, pat = a.mode, a.N, [int(x) for x in a.gens.split(",")], tuple(int(x) for x in a.pat.split(","))
    t0 = time.time()
    r = check(mode, N, gens, pat, int(a.maxstates), a.dmax, joint=a.joint)
    dead = [d for d, n in r.items() if n == 0]
    print(f"{mode} N = {N} H = <{gens}> rows {pat}: level-1 joins {r} -> "
          f"{'RULED OUT at d = ' + str(dead[0]) if dead else 'not excluded'} [{time.time() - t0:.0f}s]")


if __name__ == "__main__":
    main()

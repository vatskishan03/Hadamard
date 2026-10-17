#!/usr/bin/env python3
"""Rank the Szöllősi structures of open orders for the low-energy k-tree (search/lowe_sz.py): cyclic (Z_N) and
bordered (Z_(N-1)) structures whose exact join is out of reach (pair work > 2e10) but whose blocks can still be
enumerated (C, D pools <= MAXPOOL, A, B pools <= 1e8). Each structure/pattern with Gaussian estimate >= MINE gets the
estimate inside low-energy lists of L members (lowe_sz.run(..., estimate_only=True)); one JSON line per pattern.

    python3 search/lowe_scan.py N [N ...] [--maxpool 1.5e11] [--mine 1] [--L 2e7] [--top 8] [--out FILE]
"""
import argparse, contextlib, io, json, math, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from bnw import Struct, subgroups
import ssxx, bsz, lowe_sz


def candidates(N, bordered, maxpool, maxo=44, mino=28):
    m = N - 1 if bordered else N
    for H, gens in subgroups(m).items():
        if (m - 1) in H:
            continue
        orb = sum(math.gcd(h - 1, m) for h in H) // len(H)          # Burnside: number of H-orbits on Z_m
        if orb > maxo + 2 or orb < mino:
            continue
        st = Struct(N, gens, plain=not bordered)
        if len(st.sizesA) > maxo:
            continue
        pats = [p for _, p in bsz.patterns(st)] if bordered else ssxx.patterns(st)
        for pat in pats:
            sizes = [st.counts(k).get(r, 0) for k, r in zip(ssxx.KINDS, pat)]
            if min(sizes) == 0 or max(sizes[2:]) > maxpool or max(sizes[:2]) > 1e8:
                continue
            work = min(sizes[o[0]] * sizes[o[1]] + sizes[o[2]] * sizes[o[3]] for o in ([0, 1, 2, 3], [0, 2, 1, 3], [0, 3, 1, 2]))
            if work <= 2e10:
                continue                      # the exact sweep covers it
            yield st, gens, tuple(int(x) for x in pat), sizes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("Ns", type=int, nargs="+")
    ap.add_argument("--maxpool", type=float, default=1.5e11)
    ap.add_argument("--mine", type=float, default=1.0)
    ap.add_argument("--L", type=float, default=2e7)
    ap.add_argument("--out")
    ap.add_argument("--top", type=int, default=8)
    a = ap.parse_args()
    out = open(a.out, "a") if a.out else None
    for N in a.Ns:
        rows = []
        pre = []
        for bordered in (False, True):
            for st, gens, pat, sizes in candidates(N, bordered, a.maxpool):
                e = (bsz.gauss(st, pat, K=600) if bordered else ssxx.gauss(st, pat, K=600))
                if e == -math.inf or e < a.mine:
                    continue
                pre.append((e, bordered, st, gens, pat, sizes))
        pre.sort(key=lambda r: -r[0])
        for e, bordered, st, gens, pat, sizes in pre[:a.top]:      # only the best patterns get the list estimate
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                ein = lowe_sz.run(N, gens, pat, bordered, a.L, estimate_only=True, K=4000)
            if isinstance(ein, list):                    # ruled out before the estimate (no solution possible)
                continue
            r = {"N": N, "bordered": bordered, "gens": list(gens), "H": len(st.H), "OA": len(st.sizesA),
                 "C": st.C, "rows": list(pat), "logE": round(float(e), 2), "logEin": round(float(ein), 2),
                 "sizes": [float(s) for s in sizes]}
            rows.append(r)
            if out:
                out.write(json.dumps(r) + "\n"); out.flush()
        rows.sort(key=lambda r: -r["logEin"])
        print(N, [(r["logEin"], r["logE"], "b" if r["bordered"] else "c", r["gens"], r["rows"]) for r in rows[:3]],
              flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Production run for QH(358): the tower over R-tuples of one row-sum pattern; every hit is verified exactly and written to OUT/hits/.

For each R-tuple rank (R-tuples with distinct weight, heaviest first): stage 1 (complete A, B fibres; NS sampled C, D members;
class-level join on the Z11 columns; cached in OUT/s1_r<rank>.npz), then the (R,w)-tuples in order of their estimated fibre
product, the first TOPK of them, each resolved completely by tower.c (about 0.1 to 0.5 s).  Restartable: OUT/done.jsonl is the
checkpoint, the same command resumes.  A hit file holds the four sequences over G = Z2^4 x Z11 (+/- strings, index in mixed
radix, last factor fastest) and is the input of ../assemble.py.  Seeded: stage 1 of rank rr uses numpy default_rng(SEED * 1000 + rr).

    python3 produce.py PATTERN RANKS NS TOPK OUT SEED [--only T,T,...] [--ridx I] [--quad a,b,c,d]
      PATTERN   row sums |r_A|,|r_B|,|r_C|,|r_D|, e.g. 18,0,14,14
      RANKS     comma list of R-tuple ranks;  NS  sampled C, D members per R-fibre in stage 1 (e.g. 1e5);  TOPK  (R,w)-tuples per R-tuple
      --only    resolve only the (R,w)-tuples with these ranks (the others are skipped)
      --ridx    use the R-tuple with this index into the explicit list instead of the rank (RANKS then names only the seed offset)
      --quad    resolve only the (R,w)-tuple with these class indices (a, b, c, d) of stage 1
Needs ./pjoin and ./tower (see README).
"""
import sys, os, time, json, math, argparse
import numpy as np
from engine import *

RHO_LEAF = -13.1          # log10 hit density per leaf quadruple (rho_t 10^-24.45 divided by the Gaussian retention 10^-11.34 of lines 1,2)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pattern"); ap.add_argument("ranks"); ap.add_argument("ns"); ap.add_argument("topk"); ap.add_argument("out"); ap.add_argument("seed")
    ap.add_argument("--only", default=None); ap.add_argument("--ridx", type=int, default=None); ap.add_argument("--quad", default=None)
    a = ap.parse_args()
    pat = tuple(int(x) for x in a.pattern.split(",")); ranks = [int(x) for x in a.ranks.split(",")]
    ns = int(float(a.ns)); topk = int(a.topk); outdir = a.out; seed = int(a.seed)
    only = None if a.only is None else [int(x) for x in a.only.split(",")]
    quad_sel = None if a.quad is None else tuple(int(x) for x in a.quad.split(","))
    os.makedirs(os.path.join(outdir, "hits"), exist_ok=True); wd = os.path.join(outdir, "work"); os.makedirs(wd, exist_ok=True)
    ctx = Ctx(); S, F = ctx.S, ctx.F
    exe = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tower")
    L, tl, wts, uniq = ctx.rtuples(pat)
    print("pattern", pat, "R-tuples", len(tl), "distinct weights", len(uniq), flush=True)
    done = set(); logf = os.path.join(outdir, "done.jsonl")
    tot_cover = 0.0; tot_pred = 0.0; tot_hits = 0; ntup = 0
    if os.path.exists(logf):
        for l in open(logf):
            try:
                r_ = json.loads(l); done.add((r_["rrank"], r_["trank"])); tot_cover += r_["cover"]; tot_pred += r_["pred"]; tot_hits += r_["nhits"]; ntup += 1
            except Exception: pass
        print("resuming: %d tuples done, cover %.3g, predicted %.3g, hits %d" % (ntup, tot_cover, tot_pred, tot_hits), flush=True)
    t_start = time.time(); last = time.time()
    for rr in ranks:
        idx = uniq[rr] if a.ridx is None else a.ridx
        R = RTuple(ctx, pat, idx, L, tl)
        s1 = os.path.join(outdir, "s1_r%d.npz" % rr)
        if os.path.exists(s1):
            d = np.load(s1); quads, cl = d["quads"], [d["cl0"], d["cl1"]]
        else:
            t0 = time.time(); quads, cl = stage1(ctx, R, ns, seed * 1000 + rr, wd)
            np.savez(s1, quads=quads, cl0=cl[0], cl1=cl[1])
            print("R-tuple rank %d (idx %d, weight %.2f): stage 1 -> %d class tuples (%.0fs)" % (rr, idx, wts[idx], len(quads), time.time() - t0), flush=True)
        prf = os.path.join(outdir, "s1_prio_r%d.npy" % rr)
        if os.path.exists(prf): pr = np.load(prf)
        else: pr = tuple_priority(ctx, R, cl, quads); np.save(prf, pr)
        o = np.argsort(-pr)
        t0 = time.time(); prepare_members(ctx, R); print("   members prepared (%.0fs)" % (time.time() - t0), flush=True)
        spec = os.path.join(wd, "spec.bin"); out = os.path.join(wd, "out.txt")
        if quad_sel is not None:
            sel = [int(np.flatnonzero(o == np.flatnonzero((quads == np.array(quad_sel)).all(1))[0])[0])]
            print("   quad", quad_sel, "has (R,w)-tuple rank", sel[0], flush=True)
        else:
            sel = range(min(topk, len(o))) if only is None else [t for t in only if t < len(o)]
        for tr in sel:
            if (rr, tr) in done: continue
            t0 = time.time(); quad = quads[o[tr]]; qa, qb, qc, qd = [int(x) for x in quad]
            ia, ib = write_spec(ctx, R, spec, qa, qb, R.r[2], cl[0][qc], R.r[3], cl[1][qd])
            hits_raw, stats, l1 = run_tower(spec, out, exe)
            hits = [dict(A=R.XA[ia[h[0]]].tolist(), B=R.XB[ib[h[1]]].tolist(), C=h[2], D=h[3]) for h in hits_raw]
            cover = float(stats.split("cover")[1].split()[0]) if stats else 0.0
            pred = cover * 10 ** RHO_LEAF
            nh = 0
            for k, h in enumerate(hits):
                seqs = hit_to_sequences(ctx, h); ok = verify(ctx, seqs)
                json.dump(dict(pattern=pat, rrank=rr, ridx=idx, trank=tr, quad=[qa, qb, qc, qd], ok=bool(ok), seqs=["".join("+" if v == 1 else "-" for v in x) for x in seqs], mods=S.mods.tolist(), time=time.strftime("%Y-%m-%d %H:%M:%S")), open(os.path.join(outdir, "hits", "hit_r%d_t%d_%d.json" % (rr, tr, k)), "w"))
                print("*** HIT *** r%d t%d quad %s verified=%s" % (rr, tr, [qa, qb, qc, qd], ok), flush=True); nh += 1
            tot_cover += cover; tot_pred += pred; tot_hits += nh; ntup += 1
            rec = dict(rrank=rr, trank=tr, quad=[qa, qb, qc, qd], prio=float(pr[o[tr]]), cover=cover, pred=pred, nhits=nh, secs=round(time.time() - t0, 2))
            with open(logf, "a") as f: f.write(json.dumps(rec) + "\n")
            if time.time() - last > 120 or nh:
                last = time.time()
                print("[%.0fs] %d tuples; covered quadruples %.3g; model-predicted hits %.3g; hits found %d; last tuple %.2fs" % (time.time() - t_start, ntup, tot_cover, tot_pred, tot_hits, rec["secs"]), flush=True)
    print("done all ranks; tuples %d, cover %.3g, predicted %.3g, hits %d" % (ntup, tot_cover, tot_pred, tot_hits))

if __name__ == "__main__":
    main()

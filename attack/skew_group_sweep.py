"""Exhaustive sweep of the skew Goethals-Seidel structures over non-cyclic groups (skew_group.py) for census values.

For v with a split prime (p^2 | v, gcd(p, v/p^2) = 1), every structure of skew_group.structures (odd-order subgroups of
GL(2, p) on each split prime times an odd-order subgroup of Z_m^*) with MINO..MAXO orbits, smallest first, and every
row pattern (1, b, c, d) whose pools are not empty is searched exactly (skew_group.exact: all skew A, complete pools of
B, C and D, find4b) when its pair work |A| |largest| + |B| |C| is at most CAP, and otherwise by the k-tree of
skew_ktree.run (the densest MAXTAU split values, pools of at most MAXPOOL; it finds the 1100 solutions at the 810th split value). The Gaussian estimate is not used:
for these small structures it is unreliable in both directions (attack/NOTES.md §15, §23). Verified hits go to
$HADAMARD_WORK/skew_group_hits.jsonl in the format of `skew_group.py certify`.

  python3 attack/skew_group_sweep.py v1,v2,.. [--mino 8] [--maxo 30] [--cap 3e10] [--maxtau 2000] [--maxpool 1.2e7]
"""
import argparse, json, sys, time
from pathlib import Path as _P
sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
import skew_group as SG
import skew_ktree as K


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('vs'); ap.add_argument('--mino', type=int, default=8); ap.add_argument('--maxo', type=int, default=30)
    ap.add_argument('--cap', type=float, default=3e10)
    ap.add_argument("--maxtau", type=int, default=2000, help="k-tree fallback above the cap (0: skip)")
    ap.add_argument('--maxpool', type=float, default=1.2e7)
    a = ap.parse_args()
    log = lambda s: print(s, flush=True)
    for v in [int(x) for x in a.vs.split(',')]:
        jobs = sorted(SG.structures(v, a.mino, a.maxo), key=lambda j: j[2].O)
        log(f"=== v={v}: {len(jobs)} structures with {a.mino}..{a.maxo} orbits")
        hit = False
        for name, G, st in jobs:
            pats = SG.patterns(st, v)
            for q in pats:
                t = time.time()
                res, work = SG.exact(v, G, st, q, a.cap)
                if res is None:
                    if not a.maxtau:
                        log(f"v={v} {name} mods={list(G.mods)} O={st.O} C={st.C} rows=(1,{q}) work={work:.1e}: skipped"); continue
                    hits = K.run(v, [], maxtau=a.maxtau, top=8, only=q, sols=2, maxpool=a.maxpool, prebuilt=(G, st), verbose=False)
                    log(f"v={v} {name} mods={list(G.mods)} O={st.O} C={st.C} rows=(1,{q}) work={work:.1e}: k-tree "
                        f"({a.maxtau} values of tau, pools <= {a.maxpool:.0e}): {len(hits)} verified [{time.time() - t:.0f}s]")
                    for h in hits:
                        h.update(mods=[int(x) for x in G.mods], structure=name)
                        with open(W.hits('skew_group_hits.jsonl'), 'a') as f: f.write(json.dumps(h) + '\n')
                    if hits: hit = True; break
                    continue
                total, found = res
                log(f"v={v} {name} mods={list(G.mods)} O={st.O} C={st.C} rows=(1,{q}) work={work:.1e}: {total} solutions "
                    f"[{time.time() - t:.0f}s]")
                for had, skew, seqs in found:
                    log(f"  VERIFIED Hadamard={had} skew={skew}")
                    if had and skew:
                        rec = dict(n=v, hadamard=True, skew=True, mods=[int(x) for x in G.mods], structure=name, rows=[1, *q],
                                   seqs=[''.join('+' if x > 0 else '-' for x in s) for s in seqs])
                        with open(W.hits('skew_group_hits.jsonl'), 'a') as f: f.write(json.dumps(rec) + '\n')
                        hit = True
                if hit: break
            if hit: break
        log(f"=== v={v}: {'HIT' if hit else 'no solution in the structures searched'}")


if __name__ == '__main__':
    main()

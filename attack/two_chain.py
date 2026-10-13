"""Two-block searches (D-optimal designs, Legendre pairs) on Z_v with compressions, for composite v and larger
multiplier structures than two_circ.py can join directly.

For d | v the compressions a', b' of the two sequences to Z_d satisfy PAF_a'(s) + PAF_b'(s) = T v/d for s != 0 in Z_d
(T = 2 for D-optimal designs, -2 for Legendre pairs) and are constant on the H-orbits of Z_d. Stage 1 lists, by a
dynamic program over the orbits, every compression each block attains (with its count) and joins the two lists
exactly; stage 2 lists the members of each compatible pair of compressions (meet in the middle) and joins them on all
autocorrelation classes. Nothing is sampled; with no cap reached, a structure that yields nothing is ruled out.

  python3 attack/two_chain.py dopt|lp v gens [--d D] [--maxwork W]
  python3 attack/two_chain.py dopt|lp v --sweep [--mino 20] [--maxo 60] [--maxwork W]
Hits go to $HADAMARD_WORK/two_circ_hits.jsonl.
"""
import argparse, json, math, sys, time
from collections import defaultdict
import numpy as np
from pathlib import Path as _P
sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
import spence_lib as L
import propus as PR
import propus_chain as PC
import two_circ as TC


def search(kind, v, gens, d, a, log):
    T = 2 if kind == 'dopt' else -2
    G = L.AbGroup([v])
    st = L.Structure(G, [G.perm_from_matrix(np.array([[g % v]])) for g in gens] if gens else [G.perm_from_matrix(np.array([[1]]))])
    Hs = {1}
    while True:
        n2 = Hs | {x * g % v for x in Hs for g in gens}
        if n2 == Hs: break
        Hs = n2
    comp = PC.Comp(st, v, d, sorted({h % d for h in Hs}))
    found = []
    for ra, rb in TC.patterns(v, T):
        if min(st.pool(ra), st.pool(rb)) == 0: continue
        t0 = time.time()
        Da = comp.dp(st, ra); Db = Da if rb == ra else comp.dp(st, rb)
        if Da is None or Db is None:
            log(f"{kind} v={v} gens={gens} d={d} rows ({ra}, {rb}): compression DP above its cap; skipped"); continue
        target = np.full(d // 2, T * (v // d), np.int64)
        pb = defaultdict(list)
        for k in Db: pb[PC.cpaf(comp.full(k)).tobytes()].append(k)
        pairs = []
        for ka, na in Da.items():
            need = (target - PC.cpaf(comp.full(ka))).tobytes()
            for kb in pb.get(need, ()): pairs.append((na + Db[kb], ka, kb))
        pairs.sort(key=lambda x: -x[0])
        work = sum(p[0] for p in pairs)
        log(f"{kind} v={v} gens={gens} d={d} rows ({ra}, {rb}): orbits {st.O}, classes {st.C}; compressions {len(Da)}, "
            f"{len(Db)}; compatible pairs {len(pairs)}, members {work:.2e} (pools {st.pool(ra):.2e}, {st.pool(rb):.2e}) "
            f"[{time.time() - t0:.1f}s]")
        if work > a.maxwork:
            log(f"  skipped: members above {a.maxwork:.0e}"); continue
        if st.O > 46:
            log(f"  skipped: {st.O} orbits is too many for the meet in the middle"); continue
        tot = 0
        for i, (w_, ka, kb) in enumerate(pairs):
            YA = PC.members(st, comp, ra, ka, a.maxlist); YB = PC.members(st, comp, rb, kb, a.maxlist)
            if YA is None or YB is None: log(f"  pair {i}: list above cap; skipped"); continue
            VA = np.rint(st.class_values(YA)).astype(np.int32); VB = np.rint(st.class_values(YB)).astype(np.int32)
            sols = TC.join2(VA, VB, np.full(st.C, T, np.int32)); tot += len(sols)
            for i1, j1 in sols[:2]:
                A = YA[i1].astype(np.int64)[st.orb]; B = YB[j1].astype(np.int64)[st.orb]
                ok = TC.verify_pair(G, A, B, T)
                log(f"  hit: {'D-optimal design' if T == 2 else 'Legendre pair'} verified {ok}")
                if ok:
                    rec = dict(kind=kind, v=v, group=[v], structure=f"Z{v}<{','.join(map(str, gens))}>", gens=list(gens),
                               rows=[ra, rb], method=f'two_chain d={d}', seqs=[''.join('+' if x > 0 else '-' for x in s) for s in (A, B)])
                    with open(W.hits('two_circ_hits.jsonl'), 'a') as f: f.write(json.dumps(rec) + '\n')
                    found.append(rec)
            if found: break
        log(f"  done: {min(len(pairs), i + 1) if pairs else 0} pairs, {tot} solutions [{time.time() - t0:.0f}s]")
        if found: break
    return found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('kind', choices=['dopt', 'lp']); ap.add_argument('v', type=int); ap.add_argument('gens', nargs='?', default='')
    ap.add_argument('--d', type=int, default=0); ap.add_argument('--maxwork', type=float, default=2e9)
    ap.add_argument('--maxlist', type=float, default=3e7)
    ap.add_argument('--sweep', action='store_true'); ap.add_argument('--mino', type=int, default=20); ap.add_argument('--maxo', type=int, default=60)
    a = ap.parse_args()
    log = lambda s: print(s, flush=True)
    v = a.v; ds = [a.d] if a.d else PC.best_d(v)
    if not ds: log("v is prime: no compression"); return
    d = max(ds)
    if a.sweep:
        jobs = []
        for S, gs in PR.subgroup_gens(v).items():
            O = sum(math.gcd(h - 1, v) for h in S) // len(S)
            if a.mino <= O <= a.maxo: jobs.append((O, gs))
        jobs.sort()
        log(f"=== {a.kind} v={v}: {len(jobs)} structures with {a.mino}..{a.maxo} orbits, d={d}")
    else:
        jobs = [(0, [int(x) for x in a.gens.split(',')] if a.gens else [])]
    for _, gs in jobs:
        if search(a.kind, v, gs, d, a, log): log(f"HIT {a.kind} v={v} gens={gs}"); return
    log(f"=== {a.kind} v={v}: none")


if __name__ == '__main__':
    main()

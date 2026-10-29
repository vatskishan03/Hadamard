"""prod_5912.py -- restartable sampling production for H(5912) = 8 * 739 (README, "Order 5912"): Spence's construction over G = Z211 x F7 with blocks invariant under
H = <(19, 1)> x squares(F7*) (|H| = 45, 15 K-orbits, all row sums +1).  One round = level 1 of one allocation token (gl1.c: class-restricted join of block lists of 2000 x 200000,
at most 2000 instances), the removal of instances that are block reversals of each other, and level 2 (gl2.c) on the rest.  Round r of base seed b is a function of (b, r) alone
(rng = default_rng([b, r])): a restart repeats at most the round in progress, and repro_5912.py regenerates the round (1002, 33), the hit.
usage: prod_5912.py BASE_SEED OUTPREFIX SECS     state OUTPREFIX.state.json, one JSON line per round in OUTPREFIX.jsonl; stop at a round boundary: touch OUTPREFIX.STOP
A verified hit (K-level and exact G-level check) writes OUTPREFIX.HIT_<time>.json, is certified (gs_cert.certify_solution) and stops the run."""
import sys, os, json, time, hashlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from gs_samp import make_shape, parse_allocs, block_specs, process_instances
from gsk import row_patterns, verify_g
from gl1 import L1

SHAPE = ('spence', [211], 7, [[[19], 1]], None)      # mode, K, r, generators of H_K (units, phi), characters
TOKENS = '30.0:30.0:30.15:30.15;0.0:30.0:30.0:60.30;0.0:0.0:60.0:60.30;0.0:0.0:30.0:90.30'     # n.m of the four blocks: defect weight, m-weight (gs_samp.py)
LA, LC, MAXINST = 2000, 200000, 2000


def level1(sh, l1, allocs, pcache, rng, rnd, maxinst=MAXINST):
    """the token of round rnd and its level-1 instances (count, S, E); rng: the row pattern, the four list seeds, the hash seed"""
    (al, Rs, la, lc, ms) = allocs[rnd % len(allocs)]
    if (al, ms) not in pcache: pcache[(al, ms)] = row_patterns(sh, al, ms)
    pats = pcache[(al, ms)]
    pat = pats[int(rng.integers(0, len(pats)))] if pats else None
    seeds = [int(rng.integers(1, 1 << 62)) for _ in range(4)]
    specs = block_specs(al, Rs, la, lc, seeds, ms, sh, pattern=pat)
    specs.sort(key=lambda s: s.L)
    cnt, S, E, _ = l1.run(specs, per_class=max(30000.0, 6.0 * la), maxout=maxinst, hseed=int(rng.integers(1, 1 << 40)), stop_after=maxinst)
    return al, cnt, S, E


def dedupe(sh, S, E):
    """indices of one instance per class under the independent reversals x -> x(-.) of the four blocks (an exact symmetry of the PAF condition)"""
    key = lambda s, e: hashlib.md5(s.tobytes() + e.tobytes()).digest()[:8]
    seen, keep = set(), []
    for i in range(len(S)):
        k = tuple(min(key(S[i, b], E[i, b]), key(S[i, b][sh.negperm], E[i, b][sh.negperm])) for b in range(4))
        if k not in seen: seen.add(k); keep.append(i)
    return np.array(keep, int)


def round_instances(sh, l1, allocs, pcache, base, rnd):
    """round rnd of base seed base up to level 1: (rng for level 2, token, count, S, E of the distinct instances)"""
    rng = np.random.default_rng([base, rnd])
    al, cnt, S, E = level1(sh, l1, allocs, pcache, rng, rnd)
    keep = dedupe(sh, S, E)
    return rng, al, cnt, S[keep], E[keep]


def main():
    base = int(sys.argv[1]); outp = sys.argv[2]; secs = float(sys.argv[3])
    sh = make_shape(*SHAPE); l1 = L1(sh); allocs = parse_allocs(TOKENS, LA, LC); pcache = {}
    st_path = outp + '.state.json'
    state = json.load(open(st_path)) if os.path.exists(st_path) else dict(round=0, inst=0, uniq=0, compl=0, cpu=0.0, hits=0)
    if os.path.exists(outp + '.jsonl'):                       # a kill between the record and the state leaves a repeated round: drop it
        keep = [l for l in open(outp + '.jsonl') if l.strip() and json.loads(l)['round'] < state['round']]
        open(outp + '.jsonl', 'w').writelines(keep)
    t0 = time.time(); fo = open(outp + '.jsonl', 'a')
    print('resume at round', state['round'], flush=True)
    while time.time() - t0 < secs and not os.path.exists(outp + '.STOP'):
        r = state['round']; c0 = sum(os.times()[:4]); tw = time.time()
        rng, al, cnt, S, E = round_instances(sh, l1, allocs, pcache, base, r)
        tj = time.time() - tw
        ncomp, sols, stat, tl = process_instances(sh, S, E, rng, tl2=30.0)
        gok = []
        for bl in sols:
            ok, _ = verify_g(sh, bl)
            if ok: gok.append(bl)
            else: print('!!! G-level verification FAILED', flush=True)
        cpu = sum(os.times()[:4]) - c0
        rec = dict(t=round(time.time()), base=base, round=r, alloc=list(al), n1=int(cnt), uniq=len(S), compl=int(ncomp), verified=len(gok), t_join=round(tj, 2), t_l2=round(tl, 2), cpu=round(cpu, 1), stat=stat)
        fo.write(json.dumps(rec) + '\n'); fo.flush()
        state.update(round=r + 1, inst=state['inst'] + int(cnt), uniq=state['uniq'] + len(S), compl=state['compl'] + int(ncomp), cpu=state['cpu'] + cpu, hits=state['hits'] + len(gok))
        json.dump(state, open(st_path + '.tmp', 'w')); os.replace(st_path + '.tmp', st_path)
        print(f'[{time.time()-t0:.0f}s] r{r} {al} n1 {cnt} uniq {len(S)} compl {ncomp} verified {len(gok)} (join {tj:.1f}s l2 {tl:.1f}s) totals uniq {state["uniq"]} compl {state["compl"]} cpu {state["cpu"]:.0f}', flush=True)
        if gok:
            hit = dict(mode=sh.mode, K=sh.Kmods, r=sh.r, hg=sh.hgens_raw, T=sh.chars, alloc=list(al), round=r, base=base, blocks=[[[a.tolist() for a in blk] for blk in bl] for bl in gok[:5]])
            fn = f'{outp}.HIT_{int(time.time())}.json'
            json.dump(hit, open(fn, 'w')); print('*** HIT', fn, flush=True)
            try:
                from gs_cert import certify_solution
                print('*** certificate', certify_solution(sh, gok[0], tag=os.path.basename(outp)), flush=True)
            except Exception as ex:
                print('*** HIT but certification raised', repr(ex)[:300], flush=True)
            open(outp + '.STOP', 'w').write('hit\n'); return
    print('done', state, flush=True)


if __name__ == '__main__': main()

"""Production sampler for the REAL Goethals-Seidel-type Paley plug-ins (four free blocks): level 1 (gl1.c) + level 2 (gl2.c) + exact re-verification at G level + certification.
usage: gs_samp.py MODE Kmods r hg T allocs LA LC secs seed outprefix
   MODE  gs | spence | bord ;  Kmods json list; r odd prime power; hg json [[units, phi], ...]; T json list of characters (bord) or null
   allocs  ';'-separated list (',' also accepted when no '@' is used) of  n0[.m0]:n1[.m1]:n2[.m2]:n3[.m3][/R0:R1:R2:R3][@LA,LC]   n_i = defect weight of block i (orbit weights), m_i = weight of its m-type defects
           (s = -2e; optional but recommended: the identity fixes sum m_i = ((r+1) E - 4|K| - phi(0))/4, E = sum n_i), R_i = fixed row sum of the block
           (-1 free, 0 = row sum 0 with canonical sign, R > 0 = row sum +R, no sign gauge); blocks 0,1 are the OUTER (small, LA) lists, blocks 2,3 the INNER (large, LC) lists.
   LA, LC  list sizes (generated, regenerated from seeds, never stored); round = join of LA*LC pairs on each side.
Appends one JSON line per round to <outprefix>.jsonl; on a verified hit writes <outprefix>.HIT_<time>.json (K-level blocks, G-level sequences) and certifies.
"""
import sys, os, json, math, time, pickle
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from gsk import *
from gl1 import *
from gl2io import run_gl2


def parse_allocs(txt, LA0, LC0):
    out = []
    toks = txt.split(';') if ';' in txt else (txt.split(',') if '@' not in txt else [txt])
    for tok in toks:
        la, lc = LA0, LC0
        if '@' in tok:
            tok, ll = tok.split('@'); la, lc = [int(float(x)) for x in ll.split(',')]
        Rs = (-1, -1, -1, -1)
        if '/' in tok:
            tok, rr = tok.split('/'); Rs = tuple(int(x) for x in rr.split(':'))
        parts = tok.split(':')
        al = tuple(int(x.split('.')[0]) for x in parts)
        ms = tuple(int(x.split('.')[1]) if '.' in x else -1 for x in parts)
        assert len(al) == 4 and len(Rs) == 4
        out.append((al, Rs, la, lc, ms))
    return out


def block_specs(al, Rs, LA, LC, seeds, ms=(-1, -1, -1, -1), sh=None, pattern=None):
    """pattern: ((sigma_i, eps_i) x 4) from gsk.row_patterns: block i restricted to sum s = sigma_i, sum e = eps_i (sign gauge: sigma > 0, or sigma = 0 and eps > 0, else canonical sign)"""
    specs = []
    for b in range(4):
        R = Rs[b]
        rok = None if R < 0 else {R}
        cans = not (R > 0)
        sok = eok = None; fs_ = fe_ = None
        if pattern is not None:
            sg_, ep_ = pattern[b]
            sok = {sg_}; eok = {ep_}; fs_, fe_ = sg_, ep_
            cans = (sg_ == 0 and ep_ == 0)
        L_ = (LA if b < 2 else LC)
        if pattern is not None and sh is not None and os.environ.get('GS_EXACT_CAP'):
            capx = int(float(os.environ['GS_EXACT_CAP']))
            try:
                cntp = enum_pattern(sh, al[b], ms[b], sg_, ep_, cans, count_only=True)
                if cntp <= capx:
                    ex = enum_pattern(sh, al[b], ms[b], sg_, ep_, cans, cap=capx)
                    if ex is None: ex = (np.zeros((0, sh.N), np.int8), np.zeros((0, sh.N), np.int8))
                    specs.append(BlockSpec(al[b], S=ex[0], E=ex[1], rok=rok, cansign=cans, mw=ms[b], sok=sok, eok=eok)); continue
            except RuntimeError:
                pass
        if sh is not None and L_ * 3 > space_size(sh, al[b]) / (2 if cans else 1):
            # list would be a large fraction of the whole space: use the EXACT list instead (small spaces only)
            ex = enum_weight(sh, al[b], rok=rok, cansign=cans, mw=ms[b], sok=sok, eok=eok)
            if ex is None: ex = (np.zeros((0, sh.N), np.int8), np.zeros((0, sh.N), np.int8))
            specs.append(BlockSpec(al[b], S=ex[0], E=ex[1], rok=rok, cansign=cans, mw=ms[b], sok=sok, eok=eok)); continue
        specs.append(BlockSpec(al[b], L=(LA if b < 2 else LC), seed=int(seeds[b]), rok=rok, cansign=cans, mw=ms[b], fs=fs_, fe=fe_))
    return specs


def process_instances(sh, S, E, rng, chunk=200, tl2=30.0, maxsol=20, maxfib=1500000, log=None, stop_first=True):
    """level 2 on level-1 instances; returns (ncompl, verified solutions (K-level blocks), status counter, t_l2)"""
    t0 = time.time(); ncomp = 0; sols = []; stat = {}
    for c0 in range(0, len(S), chunk):
        sl = range(c0, min(len(S), c0 + chunk))
        insts = [[(S[i, b], E[i, b]) for b in range(4)] for i in sl]
        res, ss = run_gl2(sh, insts, maxsol=maxsol, tlimit=tl2, seed=int(rng.integers(1, 1 << 30)), maxfib=maxfib, stop_first=stop_first)
        for i, d in res.items():
            ncomp += d['compl']; stat[d['status']] = stat.get(d['status'], 0) + 1
            if d['compl'] > 0: stat['solved'] = stat.get('solved', 0) + 1
        for (i, ts) in ss:
            bl = [(insts[i][X][0].astype(np.int64), ts[X].astype(np.int64), insts[i][X][1].astype(np.int64)) for X in range(4)]
            if sh.check_klevel(bl): sols.append(bl)
            else: print('!!! level-2 solution failed the K-level check', flush=True)
    return ncomp, sols, stat, time.time() - t0


def run(sh, allocs, secs, seed, outp, log=print, certify=True, per_class=None, maxrounds=None, tl2=30.0, maxinst=100000, use_patterns=False):
    rng = np.random.default_rng(seed); l1 = L1(sh); pcache = {}; done_exact = set(); exhausted = set()
    pc_path = os.environ.get('GS_PCACHE')
    if pc_path and os.path.exists(pc_path):
        try: pcache = pickle.load(open(pc_path, 'rb'))
        except Exception: pcache = {}
    pc_new = 0
    t0 = time.time(); ai = 0; tot = dict(rounds=0, inst=0, compl=0, verified=0)
    fo = open(outp + '.jsonl', 'a')
    while time.time() - t0 < secs and (maxrounds is None or tot['rounds'] < maxrounds):
        (al, Rs, LA, LC, ms) = allocs[ai % len(allocs)]; ai += 1
        seeds = [int(rng.integers(1, 1 << 62)) for _ in range(4)]
        pat = None
        if use_patterns and all(R < 0 for R in Rs):
            key = (al, ms)
            if key not in pcache:
                tp0 = time.time(); pcache[key] = row_patterns(sh, al, ms); log(f'row patterns for {al} {ms}: {None if pcache[key] is None else len(pcache[key])} ({time.time()-tp0:.0f}s)')
                pc_new += 1
                if pc_path and pc_new % 100 == 0:
                    pickle.dump(pcache, open(pc_path + '.tmp', 'wb')); os.replace(pc_path + '.tmp', pc_path)
            if pcache[key] is not None and len(pcache[key]) == 0:
                log('no row-sum pattern for', al, ms); time.sleep(0.01); continue
            if pcache[key]: pat = pcache[key][int(rng.integers(0, len(pcache[key])))]
        specs = block_specs(al, Rs, LA, LC, seeds, ms, sh, pattern=pat)
        if all(s.S is not None for s in specs):
            # all four lists are EXACT: the round is deterministic; never repeat it
            ekey = (al, ms, Rs, pat)
            if ekey in done_exact:
                if pcache.get((al, ms)) and all(((al, ms, Rs, p_) in done_exact) for p_ in pcache[(al, ms)]) or (not pcache.get((al, ms)) and pat is None):
                    if (al, ms, Rs) not in exhausted: exhausted.add((al, ms, Rs)); log(f'allocation {al} {ms} {Rs} EXHAUSTED (exact lists, all patterns done)')
                if len(exhausted) >= len(allocs): log('all allocations exhausted'); return None, tot
                continue
            done_exact.add(ekey)
        if any(s.L == 0 for s in specs):
            log('empty list for', al, ms, Rs); ai += 0; time.sleep(0.01); continue
        specs.sort(key=lambda s: s.L)
        pc = per_class if per_class else max(30000.0, 6.0 * LA)
        tj0 = time.time(); cpu0 = sum(os.times()[:4])
        cnt, S, E, st = l1.run(specs, per_class=pc, maxout=maxinst, hseed=int(rng.integers(1, 1 << 40)), stop_after=maxinst)
        tj = time.time() - tj0; cpuj = sum(os.times()[:4]) - cpu0
        cpu1 = sum(os.times()[:4])
        ncomp, sols, stat, tl = process_instances(sh, S, E, rng, tl2=tl2)
        cpul = sum(os.times()[:4]) - cpu1
        # G-level re-verification of every solution kept
        gok = []
        for bl in sols:
            ok, seqs = verify_g(sh, bl)
            if ok: gok.append(bl)
            else: print('!!! G-level verification FAILED', flush=True)
        rec = dict(t=round(time.time()), alloc=list(al), mw=list(ms), R=list(Rs), pattern=pat, LA=LA, LC=LC, n1=int(cnt), n1_stored=len(S), compl=int(ncomp), verified=len(gok),
                   t_join=round(tj, 2), t_l2=round(tl, 2), cpu_join=round(cpuj, 2), cpu_l2=round(cpul, 2), lists=[float(x) for x in st[3:7]], pairs=[float(st[0]), float(st[1])], stat=stat, seed=seed)
        fo.write(json.dumps(rec) + '\n'); fo.flush()
        tot['rounds'] += 1; tot['inst'] += cnt; tot['compl'] += ncomp; tot['verified'] += len(gok)
        log(f'[{time.time()-t0:.0f}s] alloc {al} R {Rs} n1 {cnt} compl {ncomp} verified {len(gok)} (join {tj:.1f}s l2 {tl:.1f}s) {stat} totals inst {tot["inst"]} compl {tot["compl"]}')
        if gok:
            hit = dict(mode=sh.mode, K=sh.Kmods, r=sh.r, hg=sh.hgens_raw, T=sh.chars, alloc=list(al), R=list(Rs),
                       blocks=[[[a.tolist() for a in blk] for blk in bl] for bl in gok[:5]])
            fn = f'{outp}.HIT_{int(time.time())}.json'
            json.dump(hit, open(fn, 'w'))
            log('*** HIT', fn)
            if certify:
                try:
                    from gs_cert import certify_solution
                    log('*** certificate', certify_solution(sh, gok[0], tag=os.path.basename(outp)))
                except Exception as ex:
                    log('*** HIT but certification raised', repr(ex)[:300])
            return 0, tot
    if pc_path:
        pickle.dump(pcache, open(pc_path + '.tmp', 'wb')); os.replace(pc_path + '.tmp', pc_path)
    return None, tot


def make_shape(mode, Kmods, r, hg, T):
    hg = [(tuple(u) if isinstance(u, list) else u, ph) for u, ph in hg]
    return GShape(Kmods, r, hg, mode, [tuple(c) for c in T] if T else None)


if __name__ == '__main__':
    mode = sys.argv[1]; Kmods = json.loads(sys.argv[2]); r = int(sys.argv[3]); hg = json.loads(sys.argv[4]); T = json.loads(sys.argv[5])
    allocs = parse_allocs(sys.argv[6], int(float(sys.argv[7])), int(float(sys.argv[8])))
    secs = float(sys.argv[9]); seed = int(sys.argv[10]); outp = sys.argv[11]
    sh = make_shape(mode, Kmods, r, hg, T)
    up = os.environ.get('GS_PATTERNS', '1') == '1'
    print(run(sh, allocs, secs, seed, outp, log=lambda *a: print(*a, flush=True), use_patterns=up, maxinst=int(float(os.environ.get('GS_MAXINST', '20000')))), flush=True)

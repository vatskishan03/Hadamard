"""Pattern DFS for the shapes with H_K = 1 (every H_K-orbit is a single position): the defect sets D_X of B, C, D up to translation and reflection, A = pairs {+-a}.
Count level (valid only in the rigid regime (q-2)^2 > 2 E q + 4|K|): no lag class k with exactly one big-big pair; then types and signs by bbpat.check_pattern
(the signed bound test, always valid), then the exact row-sum pattern filter of ../sk_pat.py.  The survivors are listed; ALL patterns dead => the shape is EXACTLY empty.
usage (from search/skewplug):  python3 patterns/bbhk1.py 'Kmods' q [tlimit]      e.g.  python3 patterns/bbhk1.py "[19]" 19     (361 = 19^2: ALIVE 0)"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
import json, time, itertools, math, warnings
warnings.filterwarnings('ignore')
import numpy as np
from sk import SShape
from bbpat import PatCtx, check_pattern, identity_targets


def canon_sets(m, n, add, neg):
    """all n-subsets of K (any abelian group given by its addition table add[a][b] and negation neg) containing 0 up to translation and reflection"""
    if n == 0: return [()]
    if n == 1: return [(0,)]
    out = []
    for rest in itertools.combinations(range(1, m), n - 1):
        s = (0,) + rest
        best = None
        for t in s:
            nt = int(neg[t])
            tr = [int(add[x, nt]) for x in s]
            for refl in (False, True):
                img = tuple(sorted((int(neg[x]) if refl else x) for x in tr))
                if best is None or img < best: best = img
        if best == s: out.append(s)
    return out


def run(Kmods, q, tlimit=3600, verbose=True):
    sh = SShape(Kmods, q, []); N = sh.N
    ctx = PatCtx(sh); L = ctx.L
    tg = identity_targets(sh)
    t0 = time.time()
    # the A pairs
    pairidx = {}
    for j, (o, po) in enumerate(sh.apairs): pairidx[int(sh.spos[o][0])] = j
    # lag vectors
    def lagvec_set(D):
        d = np.zeros(N, dtype=np.int64); d[list(D)] = 1
        return (d[None, :] * d[ctx.SH]).sum(1)
    Es = sorted({ep + em for (ep, em) in tg if em is not None}) if q > 3 else []
    Etot = max(Es) if Es else 0
    rigid = (q - 2) ** 2 > 2 * Etot * q + 4 * N
    print(f'K={Kmods} q={q} N={N} targets={tg} E={Es} rigid-count-test-valid={rigid}', flush=True)
    stats = dict(leaves=0, count_alive=0, typed_tested=0, alive=[], killed_signed=0)
    cand_cache = {}
    def cands(n):
        if n not in cand_cache:
            cs = canon_sets(N, n, sh.Kg.add, sh.negperm)
            M = np.array([lagvec_set(c) for c in cs]) if cs else np.zeros((0, L), dtype=np.int64)
            cand_cache[n] = (cs, M)
        return cand_cache[n]
    def A_cands(nA):
        half = nA // 2
        if half == 0: return [()], np.zeros((1, L), dtype=np.int64)
        reps = [k for k in range(1, N) if k < int(sh.negperm[k])]; out = []; vecs = []
        for comb in itertools.combinations(reps, half):
            D = []
            for a in comb: D += [a, int(sh.negperm[a])]
            out.append(comb); vecs.append(lagvec_set(D))
        return out, np.array(vecs)
    seen_alloc = set()
    for (ep, em) in tg:
        for nA in range(0, ep + em + 1, 2):
            rem = ep + em - nA
            for nB in range(0, rem + 1):
                for nC in range(nB, rem - nB + 1):
                    nD = rem - nB - nC
                    if nD < nC: continue
                    al = (nA, nB, nC, nD)
                    if (al, ep, em) in seen_alloc: continue
                    seen_alloc.add((al, ep, em))
                    if time.time() - t0 > tlimit: stats['timeout'] = True; break
                    csA, MA = A_cands(nA)
                    csB, MB = cands(nB); csC, MC = cands(nC); csD, MD = cands(nD)
                    if not csB or not csC or not csD: continue
                    for ia, ca in enumerate(csA):
                        for ib, cb in enumerate(csB):
                            base = MA[ia] + MB[ib]
                            for ic, cc in enumerate(csC):
                                if nC == nB and ic < ib: continue
                                part = base + MC[ic]
                                tot = part[None, :] + MD
                                if rigid: ok = ~(tot == 1).any(1)
                                else: ok = np.ones(len(MD), bool)
                                idxs = np.flatnonzero(ok)
                                if nD == nC: idxs = idxs[idxs >= ic]
                                stats['leaves'] += len(MD)
                                stats['count_alive'] += len(idxs)
                                for idd in idxs:
                                    # types: all assignments of p/m to the E positions with p-weight ep, m-weight em
                                    Apos = [(a, ) for a in ca]
                                    typed = typed_alive(ctx, sh, pairidx, ca, cb, cc, csD[idd], ep, em, rigid)
                                    stats['typed_tested'] += 1
                                    if typed is None: stats['killed_signed'] += 1
                                    else: stats['alive'].append(dict(al=al, tg=(ep, em), A=ca, B=cb, C=cc, D=csD[idd], typed=typed))
                    if verbose: print(f'  alloc {al} (ep,em)={(ep,em)}: leaves {stats["leaves"]} count-alive {stats["count_alive"]} signed-alive {len(stats["alive"])} [{time.time()-t0:.0f}s]', flush=True)
    # exact row-sum pattern filter on the (allocation, m-split) classes of the survivors
    from bbpat import weights_ms
    from sk_pat import row_patterns_skew
    cache = {}; keep = []
    for a in stats['alive']:
        # convert typed bits back to (A pairs, items)
        bits = a['typed']; items = []
        for a_ in a['A']: items.append((0, a_))
        for X, cs in ((1, a['B']), (2, a['C']), (3, a['D'])):
            for x in cs: items.append((X, x))
        Apairs = []; blks = [[], [], []]
        for (X, key), t in zip(items, bits):
            if X == 0: Apairs.append((pairidx[key], t))
            else: blks[X - 1].append((int(sh.sorb[key]), t))
        al, ms = weights_ms(sh, tuple(Apairs), tuple(blks[0]), tuple(blks[1]), tuple(blks[2]))
        if (al, ms) not in cache: cache[(al, ms)] = bool(row_patterns_skew(sh, al, ms))
        if cache[(al, ms)]: keep.append(a)
    stats['alive_before_rowfilter'] = len(stats['alive']); stats['alive'] = keep
    stats['secs'] = time.time() - t0
    return stats


def typed_alive(ctx, sh, pairidx, ca, cb, cc, cd, ep, em, rigid):
    """try all type assignments (p/m per defect item; A pairs share a type) with exact weights (ep, em); return the first alive typed pattern or None"""
    items = []     # (block, key, weight)
    for a in ca: items.append((0, a, 2))
    for X, cs in ((1, cb), (2, cc), (3, cd)):
        for x in cs: items.append((X, x, 1))
    n = len(items)
    for bits in itertools.product((0, 1), repeat=n):
        wm = sum(w for (b, k, w), t in zip(items, bits) if t)
        if wm != em: continue
        wp = sum(w for (b, k, w), t in zip(items, bits) if not t)
        if wp != ep: continue
        blocks = [[] for _ in range(4)]
        for (X, key, w), t in zip(items, bits):
            if X == 0:
                j = pairidx[key]; o, po = sh.apairs[j]; blocks[0] += [(o, t, +1), (po, t, -1)]
            else: blocks[X].append((int(sh.sorb[key]), t, +1))
        ok1, ok2, wit = check_pattern(ctx, blocks, rigid=rigid)
        if ok1 and ok2: return bits
    return None


if __name__ == '__main__':
    Kmods = json.loads(sys.argv[1]); q = int(sys.argv[2]); tl = float(sys.argv[3]) if len(sys.argv) > 3 else 3600
    st = run(Kmods, q, tl)
    print('RESULT', Kmods, q, 'leaves', st['leaves'], 'count-alive', st['count_alive'], 'typed tested', st['typed_tested'], 'killed(signed)', st['killed_signed'], 'ALIVE', len(st['alive']), 'timeout', st.get('timeout', False), 'secs', round(st['secs'], 1))
    for a in st['alive'][:10]: print(a)

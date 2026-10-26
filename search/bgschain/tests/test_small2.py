import sys, os, itertools, json, tempfile, shutil, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
import small_lib as S
import hchain


def run_engine(cfg, levels, split=(), bdepth=1, sym=True, maxfibre=3_000_000, normalizer=None):
    c = dict(cfg, levels=levels, split=list(split), bdepth=bdepth, sym_swap=sym, maxfibre=maxfibre)
    if normalizer is not None: c['normalizer'] = normalizer
    d = tempfile.mkdtemp(prefix='hc_')
    try:
        ch = hchain.Chain(c, d, verbose=False)
        data = ch.prepare()
        found = []; ordered = 0; complete = True; agg = hchain.Stats()
        for b in range(len(data['w'])):
            st, ctx = ch.run_branch(data, b)
            found += ctx['found']; complete &= ctx['complete']
            for k, v in st.items(): agg.add(k, v)
            ordered += int(data['mult'][b]) * st.get('solutions_ordered_w', 0)
        return found, ordered, complete, ch, agg, len(data['w'])
    finally:
        shutil.rmtree(d)


def expand(ch, found):
    P = ch.P; out = set()
    classes = {}
    for i in range(4): classes.setdefault(tuple(P.ys[i].tolist()), []).append(i)
    # images under the normalizer
    imgs = []
    for nu in P.normalizer:
        imgs.append(np.array([P.orb[nu[P.rep_of[o]]] for o in range(P.O)]))
    allf = set(tuple(ms) for ms in found)
    frontier = list(allf)
    def act(ms, img):          # x'(g) = x(nu g):  Y'[o] = Y[img[o]]
        out = []
        for m in ms:
            b = [(int(m) >> int(img[o])) & 1 for o in range(P.O)]
            out.append(int(sum(1 << o for o in range(P.O) if b[o])))
        return tuple(out)
    while frontier:
        ms = frontier.pop()
        for img in imgs:
            m2 = act(ms, img)
            if m2 not in allf: allf.add(m2); frontier.append(m2)
    found = list(allf)
    for ms in found:
        # all assignments of the (unordered) masks to blocks with equal y
        per = []
        for key, idxs in classes.items():
            vals = [ms[i] for i in idxs]
            per.append((idxs, set(itertools.permutations(vals))))
        for combo in itertools.product(*[p[1] for p in per]):
            q = [None] * 4
            for (idxs, _), perm in zip(per, combo):
                for i, x in zip(idxs, perm): q[i] = x
            out.add(tuple(q))
    return out


def check(cfg, levels, split=(), bdepth=1, sym=True, label=''):
    sols, n = S.brute(cfg)
    bs = set(sols)
    t = time.time()
    found, ordered, complete, ch, agg, nb = run_engine(cfg, levels, split, bdepth, sym)
    es = expand(ch, found)
    ok = (es == bs) and complete and (ordered == len(bs) or not sym)
    print(f"{label:40s} brute {len(bs):6d} engine_set {len(es):6d} ordered_count {ordered:6d} branches {nb:4d} complete {complete} -> {'OK' if ok else 'MISMATCH'} [{time.time()-t:.1f}s]", flush=True)
    return ok


def check_n(cfg, levels, split=(), bdepth=1, sym=True, normalizer=None, label=''):
    sols, n = S.brute(cfg)
    bs = set(sols); t = time.time()
    found, ordered, complete, ch, agg, nb = run_engine(cfg, levels, split, bdepth, sym, normalizer=normalizer)
    es = expand(ch, found)
    ok = (es == bs) and complete and (ordered == len(bs) or not sym)
    print(f"{label:46s} brute {len(bs):6d} engine_set {len(es):6d} weighted_count {ordered:6d} branches {nb:4d} complete {complete} -> {'OK' if ok else 'MISMATCH'} [{time.time()-t:.1f}s]", flush=True)
    return ok


if __name__ == '__main__':
    allok = True
    # equal-y blocks (block symmetry) and normalizer (units of Z3 x Z5 not in H), G = Z2^2 x Z3 x Z5
    gens2 = [S.mat(S.AL, 2, 4)]          # 16 orbits
    gens1 = [S.mat(S.AL, 1, 2)]          # 12 orbits
    cases = [(gens1, 1, []), (gens1, 3, [(0, 1)]), (gens2, 21, [(0, 1)]), (gens2, 22, [(0, 1), (2, 3)]), (gens2, 24, [(0, 1), (1, 2), (2, 3)]), (gens2, 23, [(1, 2)])]
    for gens, seed, same in cases:
        cfg, seqs = S.planted([2, 2, 3, 5], 2, gens, seed, same_y=same)
        for levels, split, bd, nz in (([[3], [0, 1, 3]], [], 1, None), ([[3], [0, 1, 3]], [], 1, 'units'), ([[0, 1, 3], [0, 1, 2]], [[2, 3]], 2, 'units'), ([[2, 3]], [], 1, 'units')):
            if nz == 'units' and gens is not gens1: continue     # the planted gens2 instances are not invariant under the units normalizer (the engine refuses it at start-up)
            allok &= check_n(cfg, levels, split, bd, True, nz, label=f"seed{seed} same={same} L={levels} sp={split} nz={nz}")
        allok &= check_n(cfg, [[3], [0, 1, 3]], [], 1, False, None, label=f"seed{seed} same={same} sym=off")
    print('ALL OK' if allok else 'FAILURES')

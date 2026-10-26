#!/usr/bin/env python3
"""hchain.py -- exact compression-chain search for bordered Goethals-Seidel quadruples of any index k over an abelian group G.

    python3 hchain.py prepare CFG.json [--state DIR] [--force]       tree down to `bdepth` levels -> branch list (cached)
    python3 hchain.py run     CFG.json [--state DIR] [--worker i/n | --dynamic] [--budget S] [--stop-first] [--max-branches N]
                                       [--branches a:b] [--ids i,j] [--shuffle SEED] [--leaf-mode count|search]
    python3 hchain.py status  CFG.json [--state DIR]
    python3 hchain.py verify  CFG.json SOLUTIONS.jsonl

The state directory (default ./state/<name> under the current directory) holds branches.npz, br/<id>.done, solutions.jsonl.

Method.  Level l (levels[l-1]) = compression to a quotient of G.  Tree: the root lists, for each block, every compressed
vector on levels[0] with multiplicity (DP over the H-orbits, exact coset sums), and joins the four lists exactly on the
quotient classes; each solution is a node.  Inside a node the next level is listed by a DP restricted to the fibre of all
earlier compressions and joined again.  A branch is a node at depth `bdepth` (config), the restart unit (.done marker per
branch).  A leaf (all levels fixed) enumerates the four fibres exactly (meet in the middle over the quotient groups),
optionally splits them through the `split` levels (group by compression, join the groups first), joins the lists exactly on
all PAF classes of G, and verifies every hit over G.  Nothing is sampled; `complete` is true for a branch only when no
fibre exceeded maxfibre.
"""
import os, sys, json, time, math, argparse, hashlib
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hlib
from hlib import Problem
import modlattice


def cfg_hash(cfg):
    keys = {k: cfg[k] for k in ('mods', 'nq', 'S_gens', 'gens', 'ys', 'levels', 'bdepth', 'sym_swap', 'normalizer') if k in cfg}
    keys['tg'] = [cfg.get('t_in'), cfg.get('t_out'), cfg.get('tgt')]
    return hashlib.sha1(json.dumps(keys, sort_keys=True).encode()).hexdigest()[:12]


class StopBranch(Exception):
    pass


class Stats(dict):
    def add(self, k, v=1):
        self[k] = self.get(k, 0) + v


class Chain:
    def __init__(self, cfg, state=None, verbose=True):
        self.cfg = cfg; self.P = P = Problem(cfg)
        self.verbose = verbose
        self.name = cfg.get('name', 'run')
        self.state = state or os.path.join('state', self.name)      # default: ./state/<name> under the current directory
        os.makedirs(os.path.join(self.state, 'br'), exist_ok=True)
        self.maxfibre = int(cfg.get('maxfibre', 3_000_000))
        self.maxstates = int(cfg.get('maxstates', 60_000_000))
        self.sym = bool(cfg.get('sym_swap', True))
        self.tie_prev = [None] * 4                       # blocks with identical coset sums are interchangeable
        if self.sym:
            for i in range(1, 4):
                for j in range(i - 1, -1, -1):
                    if np.array_equal(P.ys[i], P.ys[j]): self.tie_prev[i] = j; break
        self.nlev = len(P.levels)
        # per tree level: keep only the `cap` heaviest compressed vectors per block (0 = exhaustive)
        self.caps = [int(c) for c in cfg.get('cap', [])] + [0] * self.nlev
        # per tree level: process only `child_cap` solution tuples of a node (0 = all)
        self.child_caps = [int(c) for c in cfg.get('child_cap', [])] + [0] * self.nlev
        self.bdepth = int(cfg.get('bdepth', 1))
        assert 1 <= self.bdepth <= max(self.nlev, 1)
        self.leaf_mode = cfg.get('leaf_mode', 'search')
        self.chash = cfg_hash(cfg)
        self.Ng = self._ngroup(); self._lp = {}
        self._cache = {}; self._cache_n = 0; self._cache_max = int(cfg.get('dp_cache_states', 6_000_000))
        self.offs = np.concatenate([[0], np.cumsum([P.levels[l].O for l in range(self.nlev)])]).astype(int)

    def log(self, *a):
        if self.verbose: print(*a, flush=True)

    # ------------------------------------------------------------------------------------------ tree
    def node(self, depth, fixed, tied, st, ctx, stab=None):
        """fixed: list of (Level, [4 target vectors]); tied[i]: block i is identical to its tie_prev block so far"""
        P = self.P
        if ctx.get('deadline') is not None and time.time() > ctx['deadline']:
            ctx['complete'] = False; raise StopBranch()
        if ctx['mode'] == 'collect' and depth == self.bdepth:
            vec = np.stack([np.concatenate([np.asarray(tv[i], dtype=np.int16) for _, tv in fixed]) for i in range(4)])
            ctx['collect'].append((vec, list(tied), ctx['w'], ctx.get('mult', 1)))
            return
        if depth == self.nlev:
            return self.leaf(fixed, tied, st, ctx)
        L = P.levels[depth]; d1 = depth + 1
        lists = [None] * 4; Fs = [None] * 4; t0 = time.time()
        for i in range(4):
            if tied[i]:
                lists[i] = lists[self.tie_prev[i]]; Fs[i] = Fs[self.tie_prev[i]]
            else:
                key = (depth, i, b''.join(np.asarray(tv[i], dtype=np.int16).tobytes() for _, tv in fixed))
                hit = self._cache.get(key)
                if hit is None:
                    K, cnt = P.dp(P.yorb[i], [(l, tv[i]) for l, tv in fixed], L, max_states=self.maxstates)[:2]
                    trunc = False
                    cap = self.caps[depth]
                    if cap and len(K) > cap:                       # restriction to the heaviest vectors (deterministic): NOT exhaustive
                        sel = np.sort(np.argsort(-cnt, kind='stable')[:cap]); K, cnt = K[sel], cnt[sel]; trunc = True
                    hit = (K, cnt, L.paf(K) if len(K) else np.zeros((0, 1 + L.C), dtype=np.int32), trunc)
                    self._cache[key] = hit; self._cache_n += len(K)
                    while self._cache_n > self._cache_max and len(self._cache) > 1:
                        old = next(iter(self._cache)); self._cache_n -= len(self._cache.pop(old)[0])
                else:
                    st.add('dp_cache_hits')
                if hit[3]: st.add('capped_lists'); ctx['complete'] = False
                lists[i] = hit[:2]; Fs[i] = hit[2]
        st.add('dp_time_L%d' % d1, time.time() - t0)
        st.add('dp_states_L%d' % d1, sum(len(l[0]) for l in lists))
        st.add('nodes_L%d' % d1)
        if any(len(l[0]) == 0 for l in lists):
            st.add('empty_L%d' % d1); return
        t = time.time()
        tot, sols, jst = P.join4(Fs, L.T)
        st.add('join_time_L%d' % d1, time.time() - t)
        st.add('join_pairs_L%d' % d1, float(min(len(Fs[0]) * len(Fs[1]) + len(Fs[2]) * len(Fs[3]),
                                               len(Fs[0]) * len(Fs[2]) + len(Fs[1]) * len(Fs[3]), len(Fs[0]) * len(Fs[3]) + len(Fs[1]) * len(Fs[2]))))
        if len(sols):
            sols = sols[np.lexsort((sols[:, 3], sols[:, 2], sols[:, 1], sols[:, 0]))]
            keep = np.ones(len(sols), dtype=bool)
            for i in range(4):
                j = self.tie_prev[i]
                if j is not None and tied[i]: keep &= sols[:, i] >= sols[:, j]
            sols = sols[keep]
        mult = np.ones(len(sols), dtype=np.int64)
        if stab is None: stab = list(range(len(self.Ng)))
        stabs = None
        if len(sols) and len(stab) > 1 and not self.caps[depth]:
            t = time.time(); sols, mult, stabs = self.canon(L, lists, sols, tied, stab); st.add('canon_time', time.time() - t)
        st.add('sols_L%d' % d1, len(sols))
        if len(sols) == 0: return
        w = lists[0][1][sols[:, 0]] * lists[1][1][sols[:, 1]] * lists[2][1][sols[:, 2]] * lists[3][1][sols[:, 3]]
        base_mult = ctx.get('mult', 1)
        order_r = np.argsort(-w, kind='stable')
        if self.cfg.get('child_order') == 'random' and self.child_caps[depth]:        # deterministic pseudo-random order (seed = node data), restart-safe
            seed = int.from_bytes(hashlib.sha1(b''.join(np.asarray(tv[i], dtype=np.int16).tobytes() for _, tv in fixed for i in range(4))).digest()[:8], 'little')
            order_r = np.random.default_rng(seed).permutation(len(sols))
        if self.child_caps[depth] and len(order_r) > self.child_caps[depth]:
            order_r = order_r[:self.child_caps[depth]]; ctx['complete'] = False; st.add('children_cut', len(sols) - len(order_r))
        for r in order_r:
            s = sols[r]
            ctx['mult'] = base_mult * int(mult[r])
            nf = fixed + [(L, [lists[i][0][s[i]] for i in range(4)])]
            nt = [bool(tied[i] and s[i] == s[self.tie_prev[i]]) if self.tie_prev[i] is not None else False for i in range(4)]
            ctx['w'] = float(w[r])
            self.node(depth + 1, nf, nt, st, ctx, stabs[r] if stabs is not None else stab)
            if ctx.get('stop') and ctx['found']: ctx['mult'] = base_mult; return
        ctx['mult'] = base_mult

    def _ngroup(self):
        """the group N/H generated by the normalizer, as G-permutations, identity first, one per distinct action on the H-orbits"""
        P = self.P
        key = lambda nu: tuple(int(P.orb[nu[P.rep_of[o]]]) for o in range(P.O))
        ident = np.arange(P.G.v); elems = {key(ident): ident}; frontier = [ident]
        while frontier:
            x = frontier.pop()
            for g in P.normalizer:
                y = g[x]; k = key(y)
                if k not in elems: elems[k] = y; frontier.append(y)
        ik = key(ident)
        return [elems[ik]] + [elems[k] for k in sorted(elems) if k != ik]

    def _lperms(self, L):
        if id(L) not in self._lp: self._lp[id(L)] = [self.P.level_perm(L, nu) for nu in self.Ng]
        return self._lp[id(L)]

    def stabilizer_of(self, fixed):
        """elements of N/H fixing every fixed compressed vector exactly"""
        out = []
        for gi in range(len(self.Ng)):
            ok = True
            for L, tv in fixed:
                p = self._lperms(L)[gi]
                if not all(np.array_equal(np.asarray(tv[i])[p], np.asarray(tv[i])) for i in range(4)): ok = False; break
            if ok: out.append(gi)
        return out

    def canon(self, L, lists, sols, tied, stab):
        """keep one solution tuple per orbit of the stabilizer group `stab` (indices into Ng) acting on the level coordinates (and on the
        permutations of tied blocks); mult = number of distinct tie-sorted images (orbit size); stabs = per kept solution the elements
        of stab fixing the tuple exactly (the stabilizer for the next level)."""
        P = self.P
        lp = self._lperms(L)
        classes = {}
        for i in range(4):
            j = i
            while self.tie_prev[j] is not None and tied[j]: j = self.tie_prev[j]
            classes.setdefault(j, []).append(i)
        imap = {}
        for i in range(4):
            K = lists[i][0]
            if id(K) in imap: continue
            index = {K[r].tobytes(): r for r in range(len(K))}
            imap[id(K)] = [np.array([index[np.ascontiguousarray(K[r][lp[gi]]).tobytes()] for r in range(len(K))], dtype=np.int64) for gi in stab]
        def tsort(t):
            t = list(t)
            for cl in classes.values():
                if len(cl) > 1:
                    vals = sorted(t[i] for i in cl)
                    for i, v in zip(cl, vals): t[i] = v
            return tuple(t)
        keep = []; mult = []; stabs = []
        for s in sols:
            sl = tuple(int(x) for x in s)
            imgs = []; exact_fixed = []
            for pos in range(len(stab)):
                t = tuple(int(imap[id(lists[i][0])][pos][s[i]]) for i in range(4))
                imgs.append(tsort(t))
                if t == sl: exact_fixed.append(stab[pos])
            if tsort(sl) == min(imgs):
                keep.append(s); mult.append(len(set(imgs))); stabs.append(exact_fixed)
        return (np.array(keep, dtype=np.int64).reshape(-1, 4), np.array(mult, dtype=np.int64), stabs)

    def leaf(self, fixed, tied, st, ctx):
        P = self.P; t0 = time.time()
        if ctx.get('max_leaves') is not None and ctx.setdefault('nleaf', 0) >= ctx['max_leaves']:
            ctx['complete'] = False; raise StopBranch()
        ctx['nleaf'] = ctx.get('nleaf', 0) + 1
        if self.leaf_mode == 'count':
            sizes = []
            for i in range(4):
                if tied[i]: sizes.append(sizes[self.tie_prev[i]]); continue
                sizes.append(P.fibre_count(P.yorb[i], [(l, tv[i]) for l, tv in fixed]))
            st.add('leaves'); st.add('leaves_w', ctx.get('mult', 1)); st.add('leaf_time_count', time.time() - t0)
            if min(sizes) == 0: st.add('empty_leaves'); return
            ls = sorted(math.log10(s) for s in sizes)
            st.add('leaf_sumlog10_n_i', sum(ls))
            for b, name in ((0, 'small'),):
                pass
            work = min(10 ** (ls[0] + ls[1]) + 10 ** (ls[2] + ls[3]), 10 ** (ls[0] + ls[2]) + 10 ** (ls[1] + ls[3]), 10 ** (ls[0] + ls[3]) + 10 ** (ls[1] + ls[2]))
            st.add('leaf_est_final_pairs', work)
            hist = ctx.setdefault('hist', {}); k = int(math.floor(math.log10(max(work, 1))))
            hist[k] = hist.get(k, 0) + 1
            st.add('leaf_maxn_log10_sum', ls[3])
            return
        masks = [None] * 4
        for i in range(4):
            if tied[i]:
                masks[i] = masks[self.tie_prev[i]]; continue
            fx = [(l, tv[i]) for l, tv in fixed]
            try:
                m = P.enumerate_fibre(P.yorb[i], fx)
            except OverflowError:
                st.add('truncated_leaves'); ctx['complete'] = False; return
            masks[i] = m
            if len(m) > self.maxfibre:
                st.add('truncated_leaves'); ctx['complete'] = False; return
        sizes = [len(m) for m in masks]
        st.add('leaves'); st.add('leaves_w', ctx.get('mult', 1)); st.add('leaf_list_elems', sum(sizes)); st.add('enum_time', time.time() - t0)
        ctx['t_leaf0'] = time.time()
        if min(sizes) == 0: st.add('empty_leaves'); return
        st.add('leaf_sumlog10_n_i', sum(math.log10(s) for s in sizes))
        self.final(masks, tied, st, ctx)

    def final(self, masks, tied, st, ctx):
        P = self.P; t0 = time.time()
        if P.splits:
            gid = []; gvec = []
            for i in range(4):
                if tied[i]:
                    gid.append(gid[self.tie_prev[i]]); gvec.append(gvec[self.tie_prev[i]]); continue
                u, inv = hlib.group_masks(masks[i], P.splits, P.O)
                gid.append(inv); gvec.append(u)
            Fg = []
            for i in range(4):
                parts = []; a = 0
                for sp in P.splits:
                    parts.append(sp.paf(gvec[i][:, a:a + sp.O])); a += sp.O
                Fg.append(np.concatenate(parts, axis=1))
            Tg = np.concatenate([sp.T for sp in P.splits])
            tot, gsol, _ = P.join4(Fg, Tg)
            st.add('split_time', time.time() - t0); st.add('split_groups', sum(len(g) for g in gvec)); st.add('split_sols', tot)
            quads = [s for s in gsol if not any(tied[i] and s[i] < s[self.tie_prev[i]] for i in range(4))]
            order = [np.argsort(g, kind='stable') for g in gid]
            gstart = [np.searchsorted(gid[i][order[i]], np.arange(len(gvec[i]) + 1)) for i in range(4)]
        else:
            quads = [None]
        for s in quads:
            if s is not None:
                sub = [order[i][gstart[i][s[i]]:gstart[i][s[i] + 1]] for i in range(4)]
                sm = [masks[i][sub[i]] for i in range(4)]
            else:
                sm = masks
            t = time.time()
            Vs = [P.class_values(m) for m in sm]
            st.add('final_V_time', time.time() - t); st.add('final_V_elems', sum(len(m) for m in sm))
            t = time.time()
            if os.environ.get('HCHAIN_DUMP') and ctx.setdefault('ndump', 0) < 6 and min(len(v) for v in Vs) > 1000:
                np.savez(os.path.join(os.environ['HCHAIN_DUMP'], 'dump%d.npz' % ctx['ndump']), *Vs, T=P.tgtc); ctx['ndump'] += 1
            if self.cfg.get('screen_mods'):                            # EXACT necessary condition: T - sum_b V_b[0] in the Z/p^k span of all block differences
                ts = time.time()
                passed = all(modlattice.check(Vs, P.tgtc, [tuple(m) for m in self.cfg['screen_mods']]).values())
                st.add('screen_time', time.time() - ts); st.add('screened')
                if not passed:
                    st.add('obstructed')
                    if self.cfg.get('leaflog'): ctx.setdefault('leaflog', []).append(dict(sizes=[len(v) for v in Vs], obstructed=True))
                    continue
            sc = None
            if self.cfg.get('score'):                                   # HEURISTIC ranking of the tuple (Gaussian local limit), optional pre-screen
                ts = time.time(); sc = self.score_tuple(Vs, P.tgtc); st.add('score_time', time.time() - ts); st.add('scored')
                if sc < float(self.cfg.get('score_min', -1e300)):
                    st.add('screened_out')
                    if self.cfg.get('leaflog'): ctx.setdefault('leaflog', []).append(dict(sizes=[len(v) for v in Vs], score=sc, joined=False))
                    continue
            t = time.time()
            tot, hits, stats = P.join4(Vs, P.tgtc)
            if self.cfg.get('leaflog'): ctx.setdefault('leaflog', []).append(dict(sizes=[len(v) for v in Vs], score=sc, joined=True, t_join=time.time() - t, nsol=int(tot)))
            st.add('final_join_time', time.time() - t)
            st.add('final_pairs', float(min(len(Vs[0]) * len(Vs[1]) + len(Vs[2]) * len(Vs[3]), len(Vs[0]) * len(Vs[2]) + len(Vs[1]) * len(Vs[3]),
                                            len(Vs[0]) * len(Vs[3]) + len(Vs[1]) * len(Vs[2]))))
            st.add('final_joins')
            for h in hits:
                if any(tied[i] and (s is None or s[i] == s[self.tie_prev[i]]) and h[i] < h[self.tie_prev[i]] for i in range(4)): continue
                ms = [int(sm[i][h[i]]) for i in range(4)]
                ok, msg = P.verify(ms)
                if not ok:
                    st.add('false_hits'); continue
                st.add('solutions_unordered'); st.add('solutions_ordered', self.multiplicity(ms)); st.add('solutions_ordered_w', self.multiplicity(ms) * ctx.get('mult', 1))
                ctx['found'].append(ms)
                if ctx.get('stop'): return

    @staticmethod
    def score_tuple(Vs, T):
        """Gaussian local-limit log-score of the number of solutions of V_0 + .. + V_3 = T over the four fibres (natural log, lattice covolume dropped):
        sum_b log n_b - 1/2 log det S - 1/2 d^T S^+ d, S = pooled covariance on its range, d = T - sum of the means.  -inf if d is outside the range of S
        (then the fibre sums cannot reach T).  Only a ranking heuristic: nothing about exactness depends on it."""
        mus = [v.mean(0) for v in Vs]
        S = sum(np.cov(v.T.astype(np.float64), bias=True) for v in Vs)
        d = np.asarray(T, dtype=np.float64) - sum(mus)
        w, U = np.linalg.eigh(S); keep = w > 1e-9 * w.max()
        proj = U[:, keep].T @ d
        if np.linalg.norm(d - U[:, keep] @ proj) > 1e-6 * (1.0 + np.linalg.norm(d)): return float('-inf')
        return float(sum(math.log(len(v)) for v in Vs) - 0.5 * np.log(w[keep]).sum() - 0.5 * float((proj ** 2 / w[keep]).sum()))

    def multiplicity(self, ms):
        """number of distinct orderings of the four blocks among those with equal coset sums"""
        P = self.P; groups = {}
        for i in range(4): groups.setdefault(tuple(P.ys[i].tolist()), []).append(ms[i])
        m = 1
        for g in groups.values():
            c = {}
            for x in g: c[x] = c.get(x, 0) + 1
            f = math.factorial(len(g))
            for v in c.values(): f //= math.factorial(v)
            m *= f
        return m if self.sym else 1

    # ------------------------------------------------------------------------------------------ prepare / run
    def prepare(self, force=False):
        P = self.P; path = os.path.join(self.state, 'branches.npz')
        if os.path.exists(path) and not force:
            d = np.load(path, allow_pickle=True)
            if str(d['chash']) == self.chash:
                self.log(f"branches cached: {len(d['w'])}"); return d
        ok, F, T = P.level0_check()
        self.log(f"level-0 (coset sum) check: {ok}  {F.tolist()} vs {T.tolist()}")
        if not ok: raise SystemExit("pattern is inconsistent on Q (sum of PAF of the coset sums)")
        if self.nlev == 0: raise SystemExit("no levels given")
        t0 = time.time(); st = Stats(); ctx = dict(mode='collect', collect=[], w=0.0, mult=1)
        self.node(0, [], [self.tie_prev[i] is not None for i in range(4)], st, ctx)
        col = ctx['collect']
        self.log(f"tree to depth {self.bdepth}: {len(col)} branches  " + " ".join(f"{k}={v:.4g}" for k, v in sorted(st.items())) + f" [{time.time() - t0:.1f}s]")
        if col:
            w = np.array([c[2] for c in col]); vec = np.stack([c[0] for c in col]); tied = np.array([c[1] for c in col], dtype=bool)
            mult = np.array([c[3] for c in col], dtype=np.int64)
            order = np.argsort(-w, kind='stable')
            w, vec, tied, mult = w[order], vec[order], tied[order], mult[order]
        else:
            mult = np.zeros(0, dtype=np.int64); w = np.zeros(0); vec = np.zeros((0, 4, int(self.offs[self.bdepth])), dtype=np.int16); tied = np.zeros((0, 4), dtype=bool)
        np.savez(path + '.tmp.npz', w=w, vec=vec, tied=tied, mult=mult, chash=np.array(self.chash), stats=np.array(json.dumps({k: float(v) for k, v in st.items()})))
        os.replace(path + '.tmp.npz', path)
        self.log(f"  total branch weight 10^{math.log10(w.sum()) if len(w) else float('-inf'):.2f} [{time.time() - t0:.1f}s]")
        return np.load(path, allow_pickle=True)

    def run_branch(self, d, bid, stop=False):
        P = self.P
        vec = d['vec'][bid]; tied = [bool(x) for x in d['tied'][bid]]
        fixed = []
        for l in range(self.bdepth):
            fixed.append((P.levels[l], [vec[i, self.offs[l]:self.offs[l + 1]].astype(np.int64) for i in range(4)]))
        stab = self.stabilizer_of(fixed)
        st = Stats(); ctx = dict(mode='run', found=[], complete=True, stop=stop, w=0.0, mult=1, max_leaves=self.cfg.get('max_leaves'),
                   deadline=(time.time() + self.cfg['branch_budget']) if self.cfg.get('branch_budget') else None)
        t0 = time.time()
        try:
            self.node(self.bdepth, fixed, tied, st, ctx, stab)
        except StopBranch:
            pass
        st['wall'] = time.time() - t0
        if ctx.get('leaflog'):
            os.makedirs(os.path.join(self.state, 'leaves'), exist_ok=True)
            with open(os.path.join(self.state, 'leaves', f'{bid}.jsonl'), 'w') as f:
                for r in ctx['leaflog']: f.write(json.dumps(r) + '\n')
        return st, ctx

    def sequences_json(self, ms):
        P = self.P
        seqs = ["".join('+' if x > 0 else '-' for x in P.sequence(m)) for m in ms]
        return dict(order=4 * (P.v + P.k - 1), group=[int(x) for x in P.G.mods], nq=P.nq, subgroup_generators=P.S_gens, k=P.k, gens=self.cfg['gens'], ys=[y.tolist() for y in P.ys],
                    rowsums=[int(y.sum()) for y in P.ys], autocorrelation=dict(subgroup=self.cfg.get('t_in'), elsewhere=self.cfg.get('t_out')),
                    masks=[int(m) for m in ms], sequences=seqs)

    def _claim(self, pos):
        p = os.path.join(self.state, 'br', f'{pos}.claim')
        try:
            fd = os.open(p, os.O_CREAT | os.O_EXCL | os.O_WRONLY); os.write(fd, f"{os.getpid()} {time.time()}".encode()); os.close(fd); return True
        except FileExistsError:
            return False

    def run(self, worker=None, budget=None, stop_first=False, max_branches=None, brange=None, dynamic=False, sample=None, seed=1, ids=None, shuffle=None):
        d = self.prepare()
        nb = len(d['w']); T0 = time.time(); done = 0; found = 0; agg = Stats()
        if dynamic or worker is None: positions = range(nb) if brange is None else range(brange[0], min(brange[1], nb))
        else: positions = range(worker[0], nb, worker[1]) if brange is None else range(brange[0] + worker[0], min(brange[1], nb), worker[1])
        if ids is not None: positions = ids
        if shuffle is not None:
            pset = set(positions); positions = [int(x) for x in np.random.default_rng(shuffle).permutation(nb) if int(x) in pset]
        if sample:
            rng = np.random.default_rng(seed); positions = sorted(rng.choice(nb, size=min(sample, nb), replace=False).tolist())
        for pos in positions:
            if budget is not None and time.time() - T0 > budget: break
            if max_branches is not None and done >= max_branches: break
            if stop_first and os.path.exists(os.path.join(self.state, 'FOUND')): break
            mark = os.path.join(self.state, 'br', f'{pos}.done')
            if os.path.exists(mark): continue
            if dynamic and not self._claim(pos): continue
            st, ctx = self.run_branch(d, pos, stop=stop_first)
            rec = dict(id=pos, weight=float(d['w'][pos]), orbit=int(d['mult'][pos]) if 'mult' in d.files else 1, complete=ctx['complete'], nfound=len(ctx['found']),
                       stats={k: (float(v) if isinstance(v, float) else v) for k, v in st.items()}, t=time.time())
            if ctx['found']:
                with open(os.path.join(self.state, 'solutions.jsonl'), 'a') as f:
                    for ms in ctx['found']: f.write(json.dumps(self.sequences_json(ms)) + '\n')
                found += len(ctx['found'])
                if stop_first: open(os.path.join(self.state, 'FOUND'), 'w').write(str(pos))
            with open(mark + '.tmp', 'w') as f: json.dump(rec, f)
            os.replace(mark + '.tmp', mark)
            done += 1
            for k, v in st.items(): agg.add(k, v)
            self.log(f"branch {pos}/{nb}: weight 10^{math.log10(d['w'][pos]):.1f} wall {st['wall']:.1f}s complete={ctx['complete']} found={len(ctx['found'])} "
                     + " ".join(f"{k}={v:.3g}" for k, v in sorted(st.items()) if k.startswith(('sols_L', 'leaves', 'nodes_L'))))
            if found and stop_first: break
        import resource
        self.log(f"worker: {done} branches, {found} solutions, {time.time() - T0:.1f}s, peak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.0f} MB")
        return agg, found

    def status(self):
        d = np.load(os.path.join(self.state, 'branches.npz'), allow_pickle=True)
        nb = len(d['w']); recs = []
        for f in os.listdir(os.path.join(self.state, 'br')):
            if f.endswith('.done'): recs.append(json.load(open(os.path.join(self.state, 'br', f))))
        agg = Stats()
        for r in recs:
            for k, v in r['stats'].items(): agg.add(k, v)
        print(f"{len(recs)}/{nb} branches done, {sum(r['nfound'] for r in recs)} solutions found, incomplete {sum(1 for r in recs if not r['complete'])}, "
              f"orbit-weighted ordered solution count {sum(r.get('orbit', 1) * r['stats'].get('solutions_ordered_w', r['stats'].get('solutions_ordered', 0)) for r in recs)}")
        for k, v in sorted(agg.items()): print(f"  {k}: {v:.6g}")
        # per-level rates (the root level is in prepare's stats)
        print("per-level rates over the done branches (level l = levels[l-1]):")
        for l in range(2, self.nlev + 1):
            n = agg.get('nodes_L%d' % l, 0)
            if not n: continue
            print(f"  L{l} {self.P.levels[l - 1].name}: nodes {n:.6g}, solutions/node {agg.get('sols_L%d' % l, 0) / n:.4g}, dp {1e3 * agg.get('dp_time_L%d' % l, 0) / n:.3g} ms/node, "
                  f"join {1e3 * agg.get('join_time_L%d' % l, 0) / n:.3g} ms/node, list states/node {agg.get('dp_states_L%d' % l, 0) / n:.4g}")
        lv = agg.get('leaves', 0)
        if lv:
            print(f"  leaves {lv:.6g}: mean log10(prod n_i) {agg.get('leaf_sumlog10_n_i', 0) / lv:.3g}; enum {agg.get('enum_time', 0) / lv:.3g} s/leaf, split {agg.get('split_time', 0) / lv:.3g} s/leaf, "
                  f"final V {agg.get('final_V_time', 0) / lv:.3g} s/leaf, final join {agg.get('final_join_time', 0) / lv:.3g} s/leaf")
        print(f"  wall {agg.get('wall', 0):.6g} s total over done branches")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['prepare', 'run', 'status', 'verify'])
    ap.add_argument('cfg'); ap.add_argument('extra', nargs='?')
    ap.add_argument('--state'); ap.add_argument('--worker'); ap.add_argument('--dynamic', action='store_true'); ap.add_argument('--budget', type=float)
    ap.add_argument('--stop-first', action='store_true'); ap.add_argument('--max-branches', type=int)
    ap.add_argument('--branches'); ap.add_argument('--force', action='store_true'); ap.add_argument('--leaf-mode')
    ap.add_argument('--ids'); ap.add_argument('--max-leaves', type=int); ap.add_argument('--branch-budget', type=float); ap.add_argument('--sample', type=int); ap.add_argument('--seed', type=int, default=1); ap.add_argument('--shuffle', type=int)
    a = ap.parse_args()
    cfg = json.load(open(a.cfg))
    if a.leaf_mode: cfg['leaf_mode'] = a.leaf_mode
    if a.max_leaves: cfg['max_leaves'] = a.max_leaves
    if a.branch_budget: cfg['branch_budget'] = a.branch_budget
    ch = Chain(cfg, a.state)
    if a.cmd == 'prepare': ch.prepare(force=a.force)
    elif a.cmd == 'status': ch.status()
    elif a.cmd == 'verify':
        for line in open(a.extra):
            r = json.loads(line); print(ch.P.verify(r['masks']))
    else:
        w = tuple(int(x) for x in a.worker.split('/')) if a.worker else None
        br = tuple(int(x) for x in a.branches.split(':')) if a.branches else None
        ch.run(w, a.budget, a.stop_first, a.max_branches, br, a.dynamic, a.sample, a.seed, [int(x) for x in a.ids.split(',')] if a.ids else None, a.shuffle)


if __name__ == '__main__':
    main()

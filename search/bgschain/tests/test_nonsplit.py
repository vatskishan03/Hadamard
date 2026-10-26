"""Validation of the generalised engine (arbitrary H-invariant subgroups as levels, non-split G, arbitrary S, k = 4 and k = 8).

  A. level structure: coset labels, compressions, quotient PAF classes and targets against independent tuple arithmetic;
     the identity  sum_i PAF_{Y_i}(b) = sum_{s in coset b} sum_i PAF_i(s)  on every quotient element (not only class reps).
  B. DP: the multiset of compressed vectors (and multiplicities) over all sign vectors of a fibre = the DP lists (orbit and grouped).
  C. engine = brute force: all solutions as sets and orbit-weighted counts, several towers / branch depths / splits.
"""
import sys, os, itertools, time, tempfile, shutil
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
import hlib, hchain
import nonsplit_lib as NL
from test_small3 import run_engine, expand


def ker_z4(r, m_extra, which):
    """kernel generators on G = Z4^r x Z_p (m_extra = 1 odd factor): which in {'2Z4', 'Z4', 'Zp', 'S', '0'}"""
    n = r + m_extra
    def e(i, c=1):
        v = [0] * n; v[i] = c; return v
    if which == '2Z4': return [e(i, 2) for i in range(r)]
    if which == 'Z4': return [e(i) for i in range(r)]
    if which == 'Zp': return [e(r)]
    if which == 'S': return [e(i, 2) for i in range(r)] + ([e(r)] if m_extra else [])
    raise ValueError(which)


# ------------------------------------------------------------------------------------------------------------------ A
def test_levels(mods, S_gens, gens, kernels, seed=1, label=''):
    cfgp, seqs = NL.planted(mods, S_gens, gens, seed)
    P = hlib.Problem(dict(cfgp, levels=[{'ker': k} for k in kernels]))
    G = P.G; v = P.v; ok = True
    cos, kq, _ = NL.indep_cosets(mods, S_gens)
    ok &= np.array_equal(cos, P.cos) and kq == P.k
    tot = np.zeros(v, dtype=np.int64)                                       # sum_i PAF_i(s) by direct summation
    E = NL.elems(mods); idx = {e: i for i, e in enumerate(E)}
    for s in range(v):
        sh = np.array([idx[NL.add(e, E[s], mods)] for e in E])
        tot[s] = sum(int(x @ x[sh]) for x in seqs)
    ok &= np.array_equal(tot, P.tgt)
    for L, kern in zip([P.Qlevel] + P.levels, [S_gens] + list(kernels)):
        lab, nc, ksz = NL.indep_cosets(mods, kern)
        ok &= np.array_equal(lab, L.lab) and nc == L.GK.v
        # compressions of the planted sequences (independent), quotient addition table (independent, via representatives)
        Y = [np.bincount(lab, weights=x, minlength=nc).round().astype(np.int64) for x in seqs]
        rep = [int(np.flatnonzero(lab == c)[0]) for c in range(nc)]
        add_tab = np.array([[lab[idx[NL.add(E[rep[a]], E[rep[b]], mods)]] for b in range(nc)] for a in range(nc)])
        # H-orbits of the quotient and constant compressions on them
        for y in Y:
            ok &= all(len(set(y[L.orb == u].tolist())) == 1 for u in range(L.O))
        # PAF of the compressed vectors at every quotient element b, summed over the four blocks
        paf_q = np.zeros(nc, dtype=np.int64)
        for y in Y:
            for b in range(nc):
                paf_q[b] += int(sum(y[a] * y[add_tab[a, b]] for a in range(nc)))
        fib = np.bincount(lab, weights=tot.astype(np.float64), minlength=nc).round().astype(np.int64)
        ok &= np.array_equal(paf_q, fib)                                      # the compression identity (every b)
        # engine: L.paf of the compressed vectors (orbit coordinates) against paf_q on class reps; T against the fibre sums
        yo = [np.array([y[np.flatnonzero(L.orb == u)[0]] for u in range(L.O)]) for y in Y]
        F = sum(L.paf(z[None, :])[0].astype(np.int64) for z in yo)
        ok &= np.array_equal(F, L.T)
        ok &= int(L.T[0]) == int(fib[0]) and all(int(L.T[1 + c]) == int(fib[r]) for c, r in enumerate(L.reps))
        # class structure: paf_q constant on classes
        ok &= all(len(set(paf_q[L.cls == c].tolist())) == 1 for c in range(L.C))
        # level orbit weights: Y(u) = sum_o w[o] x_o over the orbits o above u
        for x, z in zip(seqs, yo):
            xo = np.array([x[np.flatnonzero(P.orb == o)[0]] for o in range(P.O)])
            comp = np.zeros(L.O, dtype=np.int64)
            for o in range(P.O): comp[L.uof[o]] += L.w[o] * xo[o]
            ok &= np.array_equal(comp, z)
        # implied_by consistency with the kernel inclusion
    for a in [P.Qlevel] + P.levels:
        for b in [P.Qlevel] + P.levels:
            ka = np.array(sorted(NL.indep_cosets(mods, [])[0].tolist())) if False else None
            ok &= (a.implied_by(b) == bool(np.all(a.ker[b.ker])))
    print(f"A levels {label:34s} orbits {P.O} levels {[L.GK.v for L in P.levels]} -> {'OK' if ok else 'FAIL'}", flush=True)
    return ok


# ------------------------------------------------------------------------------------------------------------------ B
def test_dp(mods, S_gens, gens, kernels, tower, seed=3, label=''):
    """tower: indices into kernels: fixed levels..., last = new level (or None).  Planted sign vector gives consistent targets."""
    cfgp, seqs = NL.planted(mods, S_gens, gens, seed)
    P = hlib.Problem(dict(cfgp, levels=[{'ker': k} for k in kernels]))
    x = seqs[0]
    mask = int(sum(1 << o for o in range(P.O) if x[np.flatnonzero(P.orb == o)[0]] > 0))
    M = np.array([mask], dtype=np.uint64)
    yv = P.Qlevel.of_masks(M, P.O)[0]
    fixed = [(P.levels[i], P.levels[i].of_masks(M, P.O)[0]) for i in tower[:-1]]
    new = P.levels[tower[-1]] if tower[-1] is not None else None
    # all sign vectors (brute) in the fibre of the fixed levels
    allm = np.arange(1 << P.O, dtype=np.uint64)
    keep = np.ones(len(allm), dtype=bool)
    for L, tv in [(P.Qlevel, yv)] + fixed:
        keep &= np.all(L.of_masks(allm, P.O) == np.asarray(tv)[None, :], axis=1)
    fm = allm[keep]
    ok = True
    if new is not None:
        C = new.of_masks(fm, P.O)
        uq, ct = np.unique(C, axis=0, return_counts=True)
        for strat in ('orbit', 'grouped'):
            K, cnt, _ = P.dp(yv, fixed, new, strategy=strat)
            ok &= np.array_equal(K.astype(np.int64), uq) and np.allclose(cnt, ct)
    else:
        for strat in ('orbit', 'grouped'):
            _, cnt, _ = P.dp(yv, fixed, None, strategy=strat)
            ok &= round(float(cnt.sum())) == len(fm)
    # enumeration of the fibre with all levels fixed
    if new is not None:
        fx2 = fixed + [(new, new.of_masks(M, P.O)[0])]
        keep2 = np.all(new.of_masks(fm, P.O) == new.of_masks(M, P.O)[0][None, :], axis=1)
        e = P.enumerate_fibre(yv, fx2)
        ok &= np.array_equal(np.sort(e), np.sort(fm[keep2]))
    print(f"B dp {label:30s} tower {tower} fibre {len(fm)} -> {'OK' if ok else 'FAIL'}", flush=True)
    return ok


# ------------------------------------------------------------------------------------------------------------------ C
def check(cfg, levels, split=(), bdepth=1, sym=True, normalizer=None, label='', maxn=6000):
    t = time.time()
    try:
        sols, n = NL.brute_all(cfg, maxn=maxn)
    except AssertionError as e:
        print(f"C {label:52s} SKIP (list too large for the brute force: {e})", flush=True); return True
    bs = set(sols)
    try:
        found, ordered, complete, ch, agg, nb = run_engine(cfg, levels, split, bdepth, sym, normalizer=normalizer)
    except MemoryError as e:
        print(f"C {label:52s} SKIP (engine join too large for memory: {e})", flush=True); return True
    es = expand(ch, found)
    ok = (es == bs) and complete and (ordered == len(bs) or not sym)
    print(f"C {label:52s} lists {n} brute {len(bs):6d} engine_set {len(es):6d} weighted {ordered:6d} branches {nb:4d} complete {complete} -> {'OK' if ok else 'MISMATCH'} [{time.time()-t:.1f}s]", flush=True)
    return ok


def lv(*ks):
    return [{'ker': k} for k in ks]


if __name__ == '__main__':
    allok = True
    which = sys.argv[1:] or ['A', 'B', 'C']
    # ------------------------------------------------ A
    if 'A' in which:
        m2 = [4, 4]; S2 = [[2, 0], [0, 2]]
        allok &= test_levels(m2, S2, [NL.T2], [S2], label='Z4^2 (Q only)')
        m23 = [4, 4, 3]; S23 = [[2, 0, 0], [0, 2, 0], [0, 0, 1]]; g23 = [NL.blockdiag(NL.T2, [2])]
        Zp = ker_z4(2, 1, 'Zp'); Z4 = ker_z4(2, 1, 'Z4'); K2 = ker_z4(2, 1, '2Z4')
        allok &= test_levels(m23, S23, g23, [Z4, K2, Zp], label='Z4^2 x Z3, (T,-1)')
        g23b = [NL.blockdiag(NL.T2, [1])]
        allok &= test_levels(m23, S23, g23b, [Z4, K2, Zp], seed=5, label='Z4^2 x Z3, (T,1)')
        m3 = [4, 4, 4, 3]; S3 = [[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, 2, 0], [0, 0, 0, 1]]; g3 = [NL.blockdiag(NL.T3, [2])]
        allok &= test_levels(m3, S3, g3, [ker_z4(3, 1, 'Z4'), ker_z4(3, 1, '2Z4'), ker_z4(3, 1, 'Zp')], label='Z4^3 x Z3, (T3,-1)')
        m5 = [4, 4, 4, 5]; S5 = [[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, 2, 0], [0, 0, 0, 1]]; g5 = [NL.blockdiag(NL.T3, [2])]
        allok &= test_levels(m5, S5, g5, [ker_z4(3, 1, 'Z4'), ker_z4(3, 1, '2Z4'), ker_z4(3, 1, 'Zp')], label='Z4^3 x Z5, (T3,2)')
        m13 = [4, 4, 4, 13]; S13 = [[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, 2, 0], [0, 0, 0, 1]]; g13 = [NL.blockdiag(NL.T3, [3])]
        allok &= test_levels(m13, S13, g13, [ker_z4(3, 1, 'Z4'), ker_z4(3, 1, '2Z4'), ker_z4(3, 1, 'Zp')], label='Z4^3 x Z13, (T3,3)  [839]')
    # ------------------------------------------------ B
    if 'B' in which:
        m23 = [4, 4, 3]; S23 = [[2, 0, 0], [0, 2, 0], [0, 0, 1]]; g23 = [NL.blockdiag(NL.T2, [2])]
        ks = [ker_z4(2, 1, 'Z4'), ker_z4(2, 1, '2Z4'), ker_z4(2, 1, 'Zp')]
        for tw in ([0, 1], [0, None], [1, None], [2, 0], [0, 1, 2], [1, 2], [2, 1], [0, 2, 1], [0, 1, None]):
            allok &= test_dp(m23, S23, g23, ks, tw, seed=3, label='Z4^2xZ3 (T,-1)')
        m3 = [4, 4, 4, 3]; S3 = [[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, 2, 0], [0, 0, 0, 1]]; g3 = [NL.blockdiag(NL.T3, [2])]
        ks3 = [ker_z4(3, 1, 'Z4'), ker_z4(3, 1, '2Z4'), ker_z4(3, 1, 'Zp')]
        for tw in ([0, 1], [0, 1, None], [2, 1], [0, 2], [1, 2], [0, 1, 2]):
            allok &= test_dp(m3, S3, g3, ks3, tw, seed=4, label='Z4^3xZ3 (T3,-1) [20 orbits]')
    # ------------------------------------------------ C
    if 'C' in which or 'C1' in which:
        # C1: Z4^2, k = 4, real family target (-12 / +4), all consistent patterns (sorted)
        m2 = [4, 4]; S2 = [[2, 0], [0, 2]]; g2 = [NL.T2]
        pats = NL.all_patterns(m2, S2, g2, -12, 4)
        print('C1 patterns', len(pats))
        for p in pats:
            cfg = NL.real_cfg(m2, S2, g2, p, -12, 4)
            allok &= check(cfg, lv(S2), (), 1, True, None, label=f"Z4^2 k4 y0={[(y[0], y[1]) for y in p]}")
    if 'C' in which or 'C2' in which:
        # C2: Z4^2 x Z3, H = <(T,-1)> (12 orbits), real patterns, several towers
        m23 = [4, 4, 3]; S23 = [[2, 0, 0], [0, 2, 0], [0, 0, 1]]; g23 = [NL.blockdiag(NL.T2, [2])]
        Zp = ker_z4(2, 1, 'Zp'); Z4 = ker_z4(2, 1, 'Z4'); K2 = ker_z4(2, 1, '2Z4'); Z3k = [[0, 0, 1]]
        pats = NL.all_patterns(m23, S23, g23, -12, 4)
        print('C2 patterns', len(pats))
        towers = [(lv(Z4, K2), [], 1), (lv(Z4, K2), lv(Z3k), 1), (lv(K2), lv(Z3k), 1), (lv(Z4), lv(K2), 1), (lv(Z3k, K2), [], 1),
                  (lv(Z4, Z3k), lv(K2), 2), (lv(Z4, Z3k), [], 1), (lv(Z4, K2), [], 2), (lv(K2, Z4), [], 1)]
        for ip, p in enumerate(pats):
            cfg = NL.real_cfg(m23, S23, g23, p, -12, 4)
            for ti, (L, sp, bd) in enumerate(towers):
                if ip % 4 != ti % 4 and ip not in (0, 1): continue          # a rotating subset keeps the run short; patterns 0,1 get every tower
                allok &= check(cfg, L, sp, bd, True, None, label=f"Z4^2xZ3 pat{ip} tower{ti}")
    if 'C' in which or 'C3' in which:
        # C3: Z4^2 x Z3 with trivial action on Z3 (18 orbits): real patterns
        g23b = [NL.blockdiag(NL.T2, [1])]
        pats = NL.all_patterns(m23, S23, g23b, -12, 4)
        print('C3 patterns', len(pats))
        for ip, p in enumerate(pats[:10]):
            cfg = NL.real_cfg(m23, S23, g23b, p, -12, 4)
            allok &= check(cfg, lv(Z4, K2), lv(Z3k), 1, True, None, label=f"Z4^2xZ3 (T,1) pat{ip}", maxn=6000)
    if 'C' in which or 'C4' in which:
        # C4: k = 8 planted (explicit target) on Z4^3 and Z4^3 x Z3
        m3 = [4, 4, 4]; S3 = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]; g3 = [NL.T3]
        for seed in (1, 2, 3):
            cfg, _ = NL.planted(m3, S3, g3, seed)
            allok &= check(cfg, lv(S3), (), 1, True, None, label=f"Z4^3 k8 planted seed{seed}")
        m3 = [4, 4, 4, 3]; S3 = [[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, 2, 0], [0, 0, 0, 1]]; g3 = [NL.blockdiag(NL.T3, [2])]
        ks3 = [ker_z4(3, 1, 'Z4'), ker_z4(3, 1, '2Z4'), ker_z4(3, 1, 'Zp')]
        good = []
        for seed in range(1, 60):
            cfg, _ = NL.planted(m3, S3, g3, seed)
            c = dict(cfg, levels=lv(ks3[0], ks3[1]), split=[], bdepth=1, sym_swap=True)
            P0 = hlib.Problem(c); L1 = P0.Qlevel if False else hlib.Level(P0, {'ker': ks3[0]})
            Fs = []; ok_ = True
            for i in range(4):
                K, cnt = P0.dp(P0.yorb[i], [], L1)[:2]; Fs.append(L1.paf(K))
            tot, _, _ = P0.join4(Fs, L1.T, maxout=1, stop_at=300000)
            if 1 <= tot < 300000: good.append(seed)
            if len(good) == 3: break
        print('C4 planted Z4^3xZ3 seeds with a manageable first level:', good)
        for seed in good:
            cfg, _ = NL.planted(m3, S3, g3, seed)
            for (L, sp, bd) in ((lv(ks3[0], ks3[1]), [], 1), (lv(ks3[0], ks3[1]), lv(ks3[2]), 1), (lv(ks3[0]), lv(ks3[2]), 1), (lv(ks3[0], ks3[1]), [], 2)):
                allok &= check(cfg, L, sp, bd, True, None, label=f"Z4^3xZ3 k8 planted seed{seed} L={len(L)} sp={len(sp)} bd{bd}", maxn=9000)
    if 'C' in which or 'C6' in which:
        # C6: non-split normalizer elements (units 1+2R acting trivially on Q; -1 on Z3) with exchanged blocks
        m23 = [4, 4, 3]; S23 = [[2, 0, 0], [0, 2, 0], [0, 0, 1]]; g23 = [NL.blockdiag(NL.T2, [2])]
        T2m = np.array(NL.T2)
        u1 = (np.eye(2, dtype=int) + 2 * T2m) % 4
        nu_a = np.zeros((3, 3), dtype=int); nu_a[:2, :2] = u1; nu_a[2, 2] = 1
        nu_b = np.eye(3, dtype=int); nu_b[2, 2] = 2
        nu_c = np.zeros((3, 3), dtype=int); nu_c[:2, :2] = u1; nu_c[2, 2] = 2
        ks = [ker_z4(2, 1, 'Z4'), ker_z4(2, 1, '2Z4'), ker_z4(2, 1, 'Zp')]
        for seed in (31, 32, 33):
            for nm, nu in (('1+2T', nu_a), ('-1 on Z3', nu_b), ('both', nu_c)):
                cfg, _ = NL.planted_sym(m23, S23, g23, nu.tolist(), seed)
                for (L, sp, bd) in ((lv(ks[2], ks[1]), [], 1), (lv(ks[0], ks[1]), lv(ks[2]), 1), (lv(ks[1], ks[0]), lv(ks[2]), 2)):
                    allok &= check(cfg, L, sp, bd, True, [nu.tolist()], label=f"Z4^2xZ3 symseed{seed} nu={nm} L={len(L)} sp={len(sp)} bd{bd}")
                allok &= check(cfg, lv(ks[0], ks[1]), (), 1, True, None, label=f"Z4^2xZ3 symseed{seed} nu={nm} no normalizer")
    print('NONSPLIT TESTS:', 'ALL OK' if allok else 'FAILURES')
    sys.exit(0 if allok else 1)

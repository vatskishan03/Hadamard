"""Level-1 sampler (mitm4 on random sublists) + level-2 exact completion for orbit-constant plug-in shapes (ksym.Shape).
usage: samp.py N q b hgens(comma list u or u:phi) LC seconds seed outprefix"""
import sys, os, math, time, pickle, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ps'))
from mitm import mitm4
from pk import *
from est import sample_mode
from exact2 import enum_se


def feats_u(sh, S, E):
    out = []
    for i in range(0, len(S), 100000):
        sl = slice(i, i + 100000)
        out.append(sh.eqvec(S[sl], np.zeros_like(S[sl]), E[sl], 'u').astype(np.int32))
    return np.concatenate(out, 0)


def canon_sign(K, S, E):
    sg = np.zeros(len(S), np.int8)
    for i in range(K.nu):
        P0 = int(K.pos[i][0])
        v = np.where(S[:, P0] != 0, S[:, P0], E[:, P0])
        sg = np.where(sg == 0, v, sg)
    neg = sg < 0
    S[neg] = -S[neg]; E[neg] = -E[neg]
    return S, E


def level2(sh, blocks, TI, nU):
    """blocks: list of (s,e,Dmask,K).  returns list of solutions (each 4 x (s,t,e)) found by exact tau join"""
    from exact2 import tlists
    q = sh.q; eps = sh.eps
    pafs = np.zeros(nU, np.int64); feats = []; Tm = []
    for (s, e, Dm, K) in blocks:
        s64 = s.astype(np.int64)
        pafs += np.array([(s64 * s64[sh.shift[k]]).sum() for k in sh.repsPM])
        T_, m = tlists(sh, K, s, e, Dm)
        n_ = len(T_); S_ = np.broadcast_to(s.astype(np.int32), (n_, sh.N)); T32 = T_.astype(np.int32)
        fi = np.stack([q * (T32 * T32[:, sh.shift[k]]).sum(1) for k in sh.repsPM], 1)
        fq = np.stack([(T32 * S_[:, sh.shift[k]]).sum(1) + eps * (S_ * T32[:, sh.shift[k]]).sum(1) for k in sh.repsQ], 1) if len(sh.repsQ) else np.zeros((n_, 0), np.int32)
        feats.append(np.concatenate([fi, fq], 1).astype(np.int16)); Tm.append(T_)
    tgt = np.concatenate([TI - pafs, np.zeros(len(sh.repsQ), np.int64)])
    if np.abs(tgt).max() > 32000: return []
    n2, q2 = mitm4(feats[0], feats[1], feats[2], feats[3], tgt.astype(np.int16), logR=1, maxout=100)
    sols = []
    for (i0, i1, i2, i3) in q2[:20]:
        tt = [Tm[0][i0], Tm[1][i1], Tm[2][i2], Tm[3][i3]]
        sol = [(blocks[X][0].astype(np.int16), tt[X].astype(np.int16), blocks[X][1].astype(np.int16)) for X in range(4)]
        tot = sum(sh.eqvec(s_[None], t_[None], e_[None], 'all')[0].astype(np.int64) for s_, t_, e_ in sol)
        if (tot == sh.T.astype(np.int64)).all(): sols.append(sol)
    return sols


def main():
    N = int(sys.argv[1]); q = int(sys.argv[2]); bord = sys.argv[3] == 'b'
    hg = [tuple(int(x) for x in a.split(':')) if ':' in a else (int(a), 1) for a in sys.argv[4].split(',')] if sys.argv[4] else []
    LC = int(float(sys.argv[5])); tmax = float(sys.argv[6]); seed = int(sys.argv[7]); outp = sys.argv[8]
    allocs = [tuple(int(x) for x in a.split(':')) for a in sys.argv[9].split(',')]
    rng = np.random.default_rng(seed)
    sh = Shape(N, q, hg, bord)
    nU = sh.nU; Tu = sh.T[:nU].astype(np.int16); TI = sh.T[nU:2 * nU].astype(np.int64)
    KA, KC = sh.kAB, sh.kCD
    t0 = time.time()
    SA, EA, DA, wA = enum_se(sh, KA, True, None)
    VA = feats_u(sh, SA, EA)
    print(f'AB list {len(SA)}  [{time.time()-t0:.0f}s]', flush=True)
    stats = dict(joins=0, inst=0, filt=0, l2=0, sols=0)
    allsols = []
    ai = 0
    while time.time() - t0 < tmax:
        a = allocs[ai % len(allocs)]; ai += 1
        nA, nB, nC, nD = a
        ia = np.flatnonzero(wA == nA); ib = np.flatnonzero(wA == nB)
        if len(ia) == 0 or len(ib) == 0: continue
        def samp(n, L):
            s, t, e = sample_mode(KC, n, L, rng, 'u')
            s = s.astype(np.int8); e = e.astype(np.int8)
            s, e = canon_sign(KC, s, e)
            return s, e
        tj = time.time()
        SC, EC = samp(nC, LC); SD, ED = samp(nD, LC)
        VC = feats_u(sh, SC, EC); VD = feats_u(sh, SD, ED)
        import sympy
        R = int(sympy.nextprime(max(3, int(math.ceil((len(ia) * len(SC) + len(ib) * len(SD)) / 1.5e7)))))
        n1, quads = mitm4(VA[ia].astype(np.int16), VA[ib].astype(np.int16), VC.astype(np.int16), VD.astype(np.int16), Tu, logR=R, maxout=2_000_000, seed=int(rng.integers(1, 1 << 30)))
        stats['joins'] += 1; stats['inst'] += int(n1)
        tj2 = time.time() - tj
        # prefilter on row sums:  sum_X sigma_s^2 + q sum_X sigma_t^2 = N2 (trivial character of (I)); A,B antisymmetric t => sigma_t = 0 only if eps=-1
        nf = 0
        for (x, y, z, w) in quads:
            sC = SC[z].astype(np.int64); sD = SD[w].astype(np.int64); s_a = SA[ia[x]].astype(np.int64); s_b = SA[ib[y]].astype(np.int64)
            sig = [s_a.sum(), s_b.sum(), sC.sum(), sD.sum()]
            Rm = int(sh.N2 - sum(v * v for v in sig))
            if Rm < 0 or Rm % q: continue
            Rq = Rm // q
            # sigma_t of C,D: sum of +-1 over nondefect positions (all positions with t != 0): parity = #nondefect positions, magnitude <= that
            ndC = int((EC[z] == 0).sum()); ndD = int((ED[w] == 0).sum())
            ok = False
            for aa in range(-ndC, ndC + 1, 2):
                r2 = Rq - aa * aa
                if r2 < 0: continue
                bb = int(math.isqrt(r2))
                if bb * bb == r2 and bb <= ndD and (ndD - bb) % 2 == 0: ok = True; break
            if sh.eps == 1:
                ok = True      # A,B symmetric t: sigma_t free -> no cheap filter implemented
            if not ok: continue
            nf += 1
            stats['filt'] += 1
            # full level 2
            blocks = [(SA[ia[x]], EA[ia[x]], DA[ia[x]], KA), (SA[ib[y]], EA[ib[y]], DA[ib[y]], KA)]
            dC = np.array([bool((EC[z][K0] != 0)) for K0 in [int(p[0]) for p in KC.pos]]); dD = np.array([bool((ED[w][K0] != 0)) for K0 in [int(p[0]) for p in KC.pos]])
            blocks += [(SC[z], EC[z], dC, KC), (SD[w], ED[w], dD, KC)]
            stats['l2'] += 1
            sols = level2(sh, blocks, TI, nU)
            if sols:
                stats['sols'] += len(sols); allsols.append(sols[0])
                print('*** SOLUTION FOUND', a, flush=True)
                pickle.dump([[(x_.tolist(), y_.tolist(), z_.tolist()) for (x_, y_, z_) in s] for s in allsols], open(outp + '_sols.pkl', 'wb'))
        print(f'[{time.time()-t0:.0f}s] alloc {a}: join {tj2:.0f}s, level-1 instances {n1}, passed filter {nf}; totals {stats}', flush=True)
    print('done', stats)


if __name__ == '__main__':
    main()

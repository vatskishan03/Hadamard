"""Exhaustive two-stage solver for small-orbit plug-in shapes (ksym.Shape): stage 1 = U-system by the exact 4-list join (ps/libcmitm.so mitm4)
on (s,e,defect) configurations; stage 2 = for every level-1 instance the tau-sign system ((I) and (Q), linear hash on both) by mitm4 again.
Returns every solution of the full K-level system (up to the block-sign gauge) -> core candidates.  Used for orbCD <= ~10."""
import sys, os, math, time, itertools, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ps'))
from mitm import mitm4
from pk import *


def enum_se(sh, K, canon=True, wmax=None):
    """all (s,e,defect-units) configurations of kind K: arrays S,E (n,N) int8, defect mask (n,nu) bool; canonical block sign if canon"""
    N = sh.N
    nu = K.nu
    opts = []
    for i in range(nu):
        zu = (i == K.zero_unit)
        o = []     # (defect, s, e)
        if not K.forced[i]:
            if zu: o.append((0, 0, 0))
            else:
                o.append((0, 1, 0)); o.append((0, -1, 0))
        for be in (1, -1):
            if zu: o.append((1, -be, be))
            else:
                o.append((1, 0, be)); o.append((1, -2 * be, be))
        opts.append(o)
    sizes = [len(o) for o in opts]
    tot = int(np.prod(sizes, dtype=np.float64))
    if tot > 4e7: raise RuntimeError(f'config list too large: {tot:.3g}')
    idx = np.indices(sizes).reshape(nu, -1).T          # (tot, nu)
    S = np.zeros((idx.shape[0], N), np.int8); E = np.zeros_like(S); D = np.zeros((idx.shape[0], nu), bool)
    wt = np.zeros(idx.shape[0], np.int32)
    for i in range(nu):
        arr = np.array(opts[i], dtype=np.int8)           # (nopt, 3)
        sel = arr[idx[:, i]]
        P = K.pos[i]
        S[:, P] = sel[:, 1:2]; E[:, P] = sel[:, 2:3]; D[:, i] = sel[:, 0] == 1
        wt += (sel[:, 0] == 1) * int(K.w[i])
    keep = np.ones(len(wt), bool)
    if wmax is not None: keep &= wt <= wmax
    if canon:
        # canonical sign: first unit with nonzero s or e has s>0 or (s==0 and e>0)
        sg = np.zeros(len(wt), np.int8)
        for i in range(nu):
            P0 = int(K.pos[i][0])
            v = np.where(S[:, P0] != 0, S[:, P0], E[:, P0])
            sg = np.where(sg == 0, v, sg)
        keep &= sg >= 0
    return S[keep], E[keep], D[keep], wt[keep]


def tlists(sh, K, s, e, D):
    """tau assignments for a block with given s,e,defect mask: returns (T (n,N) int8 matrix of t-arrays, features (n, nU+nQ) int16)"""
    q = sh.q; N = sh.N
    free = [i for i in range(K.nu) if (not D[i]) and (not K.forced[i])]
    m = len(free)
    n = 1 << m
    ids = np.arange(n)
    T = np.zeros((n, N), np.int16)
    for j, i in enumerate(free):
        tau = 1 - 2 * ((ids >> j) & 1)
        P = K.pos[i]
        T[:, P] = tau[:, None] * K.sig[P][None, :]
    return T, m


def solve_tiny(sh, verbose=True, maxsol=10**9, wmax=None, logR=64, max_instances=None, maxpairs=6e10):
    t0 = time.time()
    SA, EA, DA, wA = enum_se(sh, sh.kAB, True, wmax); SC, EC, DC, wC = enum_se(sh, sh.kCD, True, wmax)
    nA, nC = len(SA), len(SC)
    if verbose: print(f'  lists A/B {nA}, C/D {nC}  (pairs {nA*nC:.3g})  [{time.time()-t0:.1f}s]', flush=True)
    def chunked(S_, E_):
        out = []
        for i in range(0, len(S_), 150000):
            sl = slice(i, i + 150000)
            out.append(sh.eqvec(S_[sl], np.zeros_like(S_[sl]), E_[sl], 'u').astype(np.int32))
        return np.concatenate(out, 0)
    VA = chunked(SA, EA); VC = chunked(SC, EC)
    if np.abs(VA).max() > 32000 or np.abs(VC).max() > 32000: raise RuntimeError('int16 overflow')
    # restrict block lists whose own k=0 energy exceeds the target
    Tu = sh.T[:sh.nU].astype(np.int64)
    emA = VA[:, 0].min(); emC = VC[:, 0].min()
    okA = VA[:, 0] <= Tu[0] - emA - 2 * emC; okC = VC[:, 0] <= Tu[0] - 2 * emA - emC
    SA, EA, DA, wA, VA = SA[okA], EA[okA], DA[okA], wA[okA], VA[okA]
    SC, EC, DC, wC, VC = SC[okC], EC[okC], DC[okC], wC[okC], VC[okC]
    nA, nC = len(SA), len(SC)
    if verbose: print(f'  after energy pruning: A/B {nA}, C/D {nC} (pairs {nA*nC:.3g})', flush=True)
    if nA * nC > maxpairs: raise RuntimeError(f'too many pairs {nA*nC:.3g}')
    if nA == 0 or nC == 0: solve_tiny.last_total_t = 0; return 0, 0, []
    logR = max(1, int(math.ceil(nA * nC / 1.5e7)))
    n1, quads = mitm4(VA.astype(np.int16), VA.astype(np.int16), VC.astype(np.int16), VC.astype(np.int16), Tu.astype(np.int16), logR=logR, maxout=20_000_000)
    if verbose: print(f'  level-1 instances (block-sign canonical): {n1}  stored {len(quads)} [{time.time()-t0:.1f}s]', flush=True)
    sols = []
    q = sh.q; eps = sh.eps; nU = sh.nU
    TI = sh.T[nU:2 * nU].astype(np.int64)
    ninst = 0; tot_t = 0
    for (a, b, c, d) in quads:
        if max_instances is not None and ninst >= max_instances: break
        ninst += 1
        blocks = [(SA[a], EA[a], DA[a], sh.kAB), (SA[b], EA[b], DA[b], sh.kAB), (SC[c], EC[c], DC[c], sh.kCD), (SC[d], EC[d], DC[d], sh.kCD)]
        # target for the t-lists:  sum q PAF_t = TI - sum PAF_s  at repsPM ; sum Q = 0
        pafs = np.zeros(nU, np.int64)
        feats = []; Tm = []
        for (s, e, Dm, K) in blocks:
            s64 = s.astype(np.int64)
            pafs += np.array([(s64 * s64[sh.shift[k]]).sum() for k in sh.repsPM])
            T_, m = tlists(sh, K, s, e, Dm)
            n_ = len(T_)
            S_ = np.broadcast_to(s.astype(np.int32), (n_, sh.N))
            T32 = T_.astype(np.int32)
            fi = np.stack([q * (T32 * T32[:, sh.shift[k]]).sum(1) for k in sh.repsPM], 1)
            fq = np.stack([(T32 * S_[:, sh.shift[k]]).sum(1) + eps * (S_ * T32[:, sh.shift[k]]).sum(1) for k in sh.repsQ], 1) if len(sh.repsQ) else np.zeros((n_, 0), np.int32)
            feats.append(np.concatenate([fi, fq], 1).astype(np.int16)); Tm.append(T_)
        tgt = np.concatenate([TI - pafs, np.zeros(len(sh.repsQ), np.int64)])
        if np.abs(tgt).max() > 32000: continue
        n2, q2 = mitm4(feats[0], feats[1], feats[2], feats[3], tgt.astype(np.int16), logR=1, maxout=1000)
        tot_t += int(n2)
        if n2 > 0:
            for (i0, i1, i2, i3) in q2[:5]:
                tt = [Tm[0][i0], Tm[1][i1], Tm[2][i2], Tm[3][i3]]
                sol = [(blocks[X][0].astype(np.int16), tt[X].astype(np.int16), blocks[X][1].astype(np.int16)) for X in range(4)]
                # full verification of the K-level system
                tot = sum(sh.eqvec(s_[None], t_[None], e_[None], 'all')[0].astype(np.int64) for s_, t_, e_ in sol)
                ok = bool((tot == sh.T.astype(np.int64)).all())
                sols.append((sol, ok))
            if verbose: print(f'  SOLUTION at instance {ninst}: {n2} t-completions (verified: {[o for _, o in sols[-5:]]})', flush=True)
            if len(sols) >= maxsol: break
    if verbose: print(f'  done: instances processed {ninst}/{n1}, solutions stored {len(sols)} [{time.time()-t0:.1f}s]', flush=True)
    solve_tiny.last_total_t = tot_t
    return n1, ninst, sols


if __name__ == '__main__':
    N = int(sys.argv[1]); q = int(sys.argv[2]); bord = sys.argv[3] == 'b'
    hg = [tuple(int(x) for x in a.split(':')) if ':' in a else (int(a), 1) for a in sys.argv[4].split(',')] if len(sys.argv) > 4 and sys.argv[4] else []
    sh = Shape(N, q, hg, bord)
    print('shape', N, q, bord, hg, 'orbCD', sh.kCD.nu, 'orbAB', sh.kAB.nu, 'R', sh.R)
    n1, ninst, sols = solve_tiny(sh)
    print('N1', n1, 'solutions stored', len(sols), 'total t-completions', solve_tiny.last_total_t)

"""Sampled two-level probe for PLAIN plug-in shapes (pk.Shape, mode 'plain'): measures the completion rate rho_a = (t-completions)/(level-1 instances)
of one defect-weight allocation, to calibrate est2 (E_a / N1_a), and returns solutions (K-level blocks) when they occur.
  level 1: A, B lists EXACT for their defect weights (canonical block sign), C, D uniform random sublists (size LC) of the weight-nC, nD (s,e)-configs;
           exact 4-list join (ps/mitm.mitm4) on the U features (PAF_u at the equation reps).
  level 2: for every level-1 instance the tau signs by an exact join on the (I) (t-part) and (Q) features.
All counts are in the canonical-sign ensemble (A, B, C, D each fixed to s-sign > 0): ordered completions per instance = same as est2's E_a/N1_a.
usage: probe.py Kmods(json) q hg(json) alloc(nA:nB:nC:nD) LC seconds seed outprefix"""
import sys, os, math, time, pickle, json
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ps'))
from mitm import mitm4
from pk import *
from est import sample_mode
from exact2 import tlists
from samp import feats_u, canon_sign


class RowFilter:
    """exact necessary condition from the trivial K-character of (I):  sum_X (sum_k s_X)^2 + q (sum_k t_X)^2 = N2.
    sum_k t_X ranges over {sum_i eps_i w_i : eps_i = +-1} over the non-defect, non-forced orbits (twist-free kinds), and is 0 for A, B when eps = -1."""
    def __init__(self, sh):
        self.sh = sh; self.cache = {}
        self.sqs = {}
    def tset(self, K, Dm, kindkey):
        key = (kindkey, Dm.tobytes())
        if key in self.cache: return self.cache[key]
        ws = [int(K.w[i]) for i in range(K.nu) if (not Dm[i]) and (not K.forced[i])]
        # twisted orbits: orbit sum of sigma may differ from the weight: use sum of sigma over the orbit
        ws = [abs(int(K.sig[K.pos[i]].sum())) for i in range(K.nu) if (not Dm[i]) and (not K.forced[i])]
        S = {0}
        for w in ws:
            S = {a + w for a in S} | {a - w for a in S}
        self.cache[key] = S
        return S
    def ok(self, sA, sB, sC, sD, DmA, DmB, DmC, DmD, KA, KC):
        sh = self.sh; q = sh.q
        Rm = int(sh.N2 - sum(int(v.sum()) ** 2 for v in (sA, sB, sC, sD)))
        if Rm < 0 or Rm % q: return False
        Rq = Rm // q
        SA = self.tset(KA, DmA, 'A'); SB = self.tset(KA, DmB, 'A')
        if sh.eps == -1: SA = {0}; SB = {0}
        SC = self.tset(KC, DmC, 'C'); SD = self.tset(KC, DmD, 'C')
        sq2 = lambda S, T: {a * a + b * b for a in S for b in T}
        P1 = sq2(SA, SB); P2 = sq2(SC, SD)
        return any((Rq - p) in P2 for p in P1)


def enum_exact_weight(sh, K, n, cap=40_000_000):
    """all (s,e) configs of kind K with defect weight exactly n, canonical sign; returns S, E (int8 (M,N)), Dmask (M,nu)"""
    nu, N = K.nu, sh.N
    w = [int(x) for x in K.w]; forced = [bool(x) for x in K.forced]
    free = [i for i in range(nu) if not forced[i]]
    base = sum(w[i] for i in range(nu) if forced[i])
    subsets = []
    def rec(j, rem, cur):
        if rem == 0:
            subsets.append(list(cur)); return
        if j == len(free) or rem < 0: return
        i = free[j]
        if w[i] <= rem:
            cur.append(i); rec(j + 1, rem - w[i], cur); cur.pop()
        rec(j + 1, rem, cur)
    if n - base < 0: return None
    rec(0, n - base, [])
    NOND = np.array([[1, 0], [-1, 0]], dtype=np.int8)
    DEF = np.array([[0, 1], [-2, 1], [0, -1], [2, -1]], dtype=np.int8)          # (s, e): p: s = 0 ; m: s = -2e
    outS, outE, outD = [], [], []
    tot = 0
    for sub in subsets:
        dset = set(sub) | {i for i in range(nu) if forced[i]}
        opts = [DEF if i in dset else NOND for i in range(nu)]
        sizes = [len(o) for o in opts]
        M = int(np.prod(sizes, dtype=np.float64)); tot += M
        if tot > cap: raise RuntimeError('AB list too large')
        idx = np.indices(sizes).reshape(nu, -1).T
        S = np.zeros((M, N), np.int8); E = np.zeros((M, N), np.int8); D = np.zeros((M, nu), bool)
        for i in range(nu):
            sel = opts[i][idx[:, i]]
            P = K.pos[i]
            S[:, P] = sel[:, 0:1]; E[:, P] = sel[:, 1:2]
            D[:, i] = (i in dset)
        # canonical sign
        sg = np.zeros(M, np.int8)
        for i in range(nu):
            P0 = int(K.pos[i][0])
            v = np.where(S[:, P0] != 0, S[:, P0], E[:, P0])
            sg = np.where(sg == 0, v, sg)
        keep = sg > 0
        outS.append(S[keep]); outE.append(E[keep]); outD.append(D[keep])
    if not outS: return None
    return np.concatenate(outS), np.concatenate(outE), np.concatenate(outD)


def level2_count(sh, blocks, TI, nU, want=2):
    """exact number of tau completions of a level-1 instance (blocks = (s,e,Dmask,K) x4); returns (count, solutions)"""
    q = sh.q; eps = sh.eps
    pafs = np.zeros(nU, np.int64); feats = []; Tm = []; sizes = []
    for (s, e, Dm, K) in blocks:
        s64 = s.astype(np.int64)
        pafs += np.array([(s64 * s64[sh.shift[k]]).sum() for k in sh.repsPM])
        T_, m = tlists(sh, K, s, e, Dm)
        n_ = len(T_); S_ = np.broadcast_to(s.astype(np.int32), (n_, sh.N)); T32 = T_.astype(np.int32)
        fi = np.stack([q * (T32 * T32[:, sh.shift[k]]).sum(1) for k in sh.repsPM], 1)
        fq = np.stack([(T32 * S_[:, sh.shift[k]]).sum(1) + eps * (S_ * T32[:, sh.shift[k]]).sum(1) for k in sh.repsQ], 1) if len(sh.repsQ) else np.zeros((n_, 0), np.int32)
        feats.append(np.concatenate([fi, fq], 1)); Tm.append(T_); sizes.append(n_)
    tgt = np.concatenate([TI - pafs, np.zeros(len(sh.repsQ), np.int64)])
    if np.abs(tgt).max() > 30000 or any(np.abs(f).max() > 30000 for f in feats): return 0, []
    R = int(max(1, math.ceil(sizes[0] * sizes[2] / 1.5e7)))
    n2, q2 = mitm4(feats[0].astype(np.int16), feats[1].astype(np.int16), feats[2].astype(np.int16), feats[3].astype(np.int16), tgt.astype(np.int16), logR=R, maxout=50)
    sols = []
    for (i0, i1, i2, i3) in q2[:want]:
        tt = [Tm[0][i0], Tm[1][i1], Tm[2][i2], Tm[3][i3]]
        sol = [(blocks[X][0].astype(np.int16), tt[X].astype(np.int16), blocks[X][1].astype(np.int16)) for X in range(4)]
        tot = sum(sh.eqvec(s_[None], t_[None], e_[None], 'all')[0].astype(np.int64) for s_, t_, e_ in sol)
        if (tot == sh.T.astype(np.int64)).all(): sols.append(sol)
    return int(n2), sols


def main():
    Kmods = json.loads(sys.argv[1]); q = int(sys.argv[2])
    hg = [(tuple(u) if isinstance(u, list) else u, ph) for u, ph in json.loads(sys.argv[3])]
    alloc = tuple(int(x) for x in sys.argv[4].split(':')); LC = int(float(sys.argv[5])); tmax = float(sys.argv[6]); seed = int(sys.argv[7]); outp = sys.argv[8]
    rng = np.random.default_rng(seed)
    sh = Shape(Kmods, q, hg, mode='plain')
    nU = sh.nU; Tu = sh.T[:nU].astype(np.int16); TI = sh.T[nU:2 * nU].astype(np.int64)
    KA, KC = sh.kAB, sh.kCD
    nA, nB, nC, nD = alloc
    t0 = time.time()
    LA = enum_exact_weight(sh, KA, nA); LB = LA if nB == nA else enum_exact_weight(sh, KA, nB)
    if LA is None or LB is None: print('empty AB list'); return
    VA = feats_u(sh, LA[0], LA[1]); VB = VA if LB is LA else feats_u(sh, LB[0], LB[1])
    print(f'shape K={Kmods} q={q} H={hg} alloc {alloc}: orbAB {KA.nu} orbCD {KC.nu}; AB lists {len(VA)}, {len(VB)} [{time.time()-t0:.0f}s]', flush=True)
    tot = dict(rounds=0, inst=0, complet=0, solvable=0, tl2=0.0, tj=0.0)
    RF = RowFilter(sh); USEFILT = (os.environ.get('NOFILT') is None)
    sols_all = []
    import sympy
    while time.time() - t0 < tmax:
        tj = time.time()
        def samp(n, L):
            s, t, e = sample_mode(KC, n, L, rng, 'u'); s = s.astype(np.int8); e = e.astype(np.int8)
            return canon_sign(KC, s, e)
        SC, EC = samp(nC, LC); SD, ED = samp(nD, LC)
        VC = feats_u(sh, SC, EC); VD = feats_u(sh, SD, ED)
        R = int(sympy.nextprime(max(3, int(math.ceil((len(VA) * len(SC) + len(VB) * len(SD)) / 1.5e7)))))
        n1, quads = mitm4(VA.astype(np.int16), VB.astype(np.int16), VC.astype(np.int16), VD.astype(np.int16), Tu, logR=R, maxout=200000, seed=int(rng.integers(1, 1 << 30)))
        tot['tj'] += time.time() - tj
        tl = time.time(); ncomp = 0; nsolv = 0; ninst = 0
        for (x, y, z, w_) in quads:
            if time.time() - t0 > tmax: break
            ninst += 1
            def dm(E_, K):
                return np.array([bool(E_[int(p[0])] != 0) for p in K.pos])
            dC = dm(EC[z], KC); dD = dm(ED[w_], KC)
            passed = RF.ok(LA[0][x], LB[0][y], SC[z], SD[w_], LA[2][x], LB[2][y], ~dC, ~dD, KA, KC) if False else RF.ok(LA[0][x], LB[0][y], SC[z], SD[w_], LA[2][x], LB[2][y], dC, dD, KA, KC)
            tot['filt_pass'] = tot.get('filt_pass', 0) + int(passed)
            if USEFILT and not passed: continue
            tot['l2_inst'] = tot.get('l2_inst', 0) + 1
            blocks = [(LA[0][x], LA[1][x], LA[2][x], KA), (LB[0][y], LB[1][y], LB[2][y], KA), (SC[z], EC[z], dC, KC), (SD[w_], ED[w_], dD, KC)]
            c, sols = level2_count(sh, blocks, TI, nU)
            if c > 0 and not passed: tot['filter_violations'] = tot.get('filter_violations', 0) + 1
            ncomp += c; nsolv += (c > 0)
            if sols:
                sols_all.append(sols[0])
                pickle.dump([[(a.tolist(), b.tolist(), c_.tolist()) for (a, b, c_) in s] for s in sols_all], open(outp + '_sols.pkl', 'wb'))
                print('*** SOLUTION', alloc, flush=True)
        tot['rounds'] += 1; tot['inst'] += ninst; tot['complet'] += ncomp; tot['solvable'] += nsolv; tot['tl2'] += time.time() - tl
        print(f'[{time.time()-t0:.0f}s] join {time.time()-tj:.1f}s level-1 instances in join {n1} processed {ninst}; totals inst {tot["inst"]} filt_pass {tot.get("filt_pass",0)} l2_inst {tot.get("l2_inst",0)} completions {tot["complet"]} solvable {tot["solvable"]} viol {tot.get("filter_violations",0)} (L2 {tot["tl2"]:.0f}s, join {tot["tj"]:.0f}s)', flush=True)
    print('done', tot, 'rho_measured(completions/inst) =', tot['complet'] / max(1, tot['inst']))


if __name__ == '__main__':
    main()

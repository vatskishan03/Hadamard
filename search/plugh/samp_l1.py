"""Production sampler (C level 1 + C level 2) for PLAIN Paley plug-in shapes with cyclic K.
 level 1 (l1.c): exact A,B lists x LC random C,D configs per round (regenerated from seeds), class-sorted hash join on the U features;
 level 2 (pl2.c): exact t-sign completion of every level-1 instance (row-sum filter, PSD, MITM).
Every completion is verified in python against the K-level equations (pk.Shape.eqvec == T) and certified at once (certify.py: G-level sequences, K K* = 2|G| I,
verify_complex.py of the repository given by HADAMARD_REPO).
usage: samp_l1.py Kmods(json) q hg(json) allocs(nA:nB:nC:nD,...) LC secs seed outprefix
Appends one JSON line per round to <outprefix>.jsonl; on a verified hit writes <outprefix>.FOUND (json with the certificate info) and <outprefix>_sols.pkl and exits 0."""
import sys, os, json, math, time, pickle
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, 'ps'))
from pk import *
from probe import enum_exact_weight
from l1 import KindC
from l2io import run_pl2


def run(Kmods, q, hg, allocs, LC, secs, seed, outp, chunk=1500, log=print):
    rng = np.random.default_rng(seed)
    sh = Shape(Kmods, q, hg, mode='plain'); kc = KindC(sh)
    KA = sh.kAB
    cacheA = {}
    def getA(n):
        if n not in cacheA: cacheA[n] = enum_exact_weight(sh, KA, n)
        return cacheA[n]
    cacheR = {}
    def getAred(n):
        # A reduced modulo the automorphisms of K (all blocks move together: every solution orbit has a representative with A in the reduced list); only for H = 1
        if n not in cacheR:
            L = getA(n)
            if L is None or len(sh.hgens) > 0 or len(sh.Kmods) != 1: cacheR[n] = L
            else:
                from autred import reduce_list
                S_, E_ = reduce_list(sh, KA, L[0], L[1])
                cacheR[n] = (S_, E_, None)
        return cacheR[n]
    t0 = time.time(); ai = 0; tot = dict(rounds=0, inst=0, compl=0)
    fo = open(outp + '.jsonl', 'a')
    while time.time() - t0 < secs:
        a = allocs[ai % len(allocs)]; ai += 1
        LCa = LC
        if isinstance(a[-1], float): LCa = int(a[-1]); a = a[:-1]
        nA, nB, nC, nD = a
        LA = getAred(nA); LB = getA(nB)
        if LA is None or LB is None: log('empty A/B list for', a); continue
        tj0 = time.time()
        rs = int(rng.integers(1, 1 << 40))
        # balance the two join sides: |A_red| * LC_eff = |B| * LD_eff with LC_eff * LD_eff = LCa^2
        ratio = max(1.0, len(LB[0]) / max(1, len(LA[0])))
        LCe = int(LCa * math.sqrt(ratio)); LDe = int(LCa / math.sqrt(ratio))
        cnt, sC, eC, sD, eD, ab = kc.run(nC, nD, LCe, rs, LA[0], LA[1], LB[0], LB[1], LD=LDe)
        tj = time.time() - tj0; tl0 = time.time()
        ncomp = 0; sols = []
        for c0 in range(0, len(ab), chunk):
            sl = slice(c0, c0 + chunk)
            insts = [[(LA[0][x], LA[1][x]), (LB[0][y], LB[1][y]), (sC[i], eC[i]), (sD[i], eD[i])] for i, (x, y) in zip(range(c0, c0 + len(ab[sl])), ab[sl])]
            res, ss = run_pl2(sh, insts, maxsol=20, seed=int(rng.integers(1, 1 << 30)))
            for i, d in res.items(): ncomp += d['compl']
            for (i, ts_) in ss:
                sol = [(insts[i][X][0].astype(np.int16), ts_[X], insts[i][X][1].astype(np.int16)) for X in range(4)]
                tot_eq = sum(sh.eqvec(s_[None], t_[None], e_[None], 'all')[0].astype(np.int64) for s_, t_, e_ in sol)
                if (tot_eq == sh.T.astype(np.int64)).all(): sols.append(sol)
        tl = time.time() - tl0
        rec = dict(t=round(time.time()), alloc=list(a), LC=LCa, LCe=LCe, LDe=LDe, n1=int(cnt), compl=int(ncomp), verified=len(sols), t_join=round(tj, 1), t_l2=round(tl, 1), seed=seed, rs=rs)
        fo.write(json.dumps(rec) + '\n'); fo.flush()
        tot['rounds'] += 1; tot['inst'] += cnt; tot['compl'] += ncomp
        log(f'[{time.time()-t0:.0f}s] alloc {a} n1 {cnt} compl {ncomp} verified {len(sols)} (join {tj:.0f}s l2 {tl:.0f}s) totals inst {tot["inst"]} compl {tot["compl"]}')
        if sols:
            pk_ = outp + '_sols_%d.pkl' % int(time.time())
            pickle.dump([[(x.tolist(), y.tolist(), z.tolist()) for (x, y, z) in s] for s in sols], open(pk_, 'wb'))
            cert = None
            try:
                from certify import certify
                cert = certify(sh, sols[0], tag=str(time.time()))
                log('*** HIT', cert)
            except Exception as ex:
                cert = dict(error=repr(ex)[:300]); log('*** HIT but certification raised', cert)
            open(outp + '.FOUND', 'w').write(json.dumps(dict(cert=cert, alloc=list(a), K=Kmods, q=q, hg=hg, pickle=pk_, n_solutions=len(sols)), default=str))
            return 0
    return tot


if __name__ == '__main__':
    Kmods = json.loads(sys.argv[1]); q = int(sys.argv[2])
    hg = [(tuple(u) if isinstance(u, list) else u, ph) for u, ph in json.loads(sys.argv[3])]
    allocs = []
    for tok in sys.argv[4].split(','):
        lc = None
        if '@' in tok: tok, lc = tok.split('@'); lc = float(lc)
        al = tuple(int(x) for x in tok.split(':'))
        allocs.append(al + ((lc,) if lc else ()))
    LC = int(float(sys.argv[5])); secs = float(sys.argv[6]); seed = int(sys.argv[7]); outp = sys.argv[8]
    print(run(Kmods, q, hg, allocs, LC, secs, seed, outp, log=lambda *a: print(*a, flush=True)), flush=True)

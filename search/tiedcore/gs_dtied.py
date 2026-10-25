"""Doubly tied Goethals-Seidel census: (X, X o n^-1, Y, Y o n^-1), X, Y invariant under H, sum PAF = 0 off 0.
Usage: gs_dtied.py maxC log spec ..."""
import sys, os, time, json, numpy as np
from gengroup import *

def gs_dtied(G, Hgens, maxC, stop_at=1):
    m = G.m
    sigs = [G.perm(h) for h in Hgens]
    OC, _ = G.orbits(sigs)
    bC = len(OC)
    if bC > maxC: return None, bC
    eng = Engine(G, sigs)
    XC = all_signs(bC); PC = eng.paf(XC, OC)
    K = len(eng.key_reps); R = keyvec(K)
    target = np.zeros(K, dtype=np.int64); target[eng.labk[0]] = 4 * m
    kT = int(target @ R)
    H = G.closure(Hgens)
    seen = set(); total = 0; first = None; ncls = 0
    for n in G.mult_elems():
        cls = frozenset(G.mcomp(n, h) for h in H)
        if cls in seen or n in H: continue
        seen.add(cls); ncls += 1
        ninv = tuple(f.minv(x) for f, x in zip(G.f, n))
        kp = eng.key_perm(G.perm(ninv))
        Rd = np.zeros_like(R); np.add.at(Rd, kp, R)
        kS = lkeys(PC, R + Rd)
        for i, j in match(kS, kT - kS):
            if i > j: continue
            if np.array_equal(PC[i].astype(np.int64) + PC[i][kp] + PC[j] + PC[j][kp], target):
                total += 1
                if first is None:
                    X = XC[i] @ OC; Y = XC[j] @ OC; P = G.perm(n)
                    X2 = np.zeros(m, dtype=np.int64); X2[P] = X; Y2 = np.zeros(m, dtype=np.int64); Y2[P] = Y
                    shape = tuple(s for f in G.f for s in f.shape); tot = 0
                    for f in (X, X2, Y, Y2):
                        Fh = np.fft.fftn(np.asarray(f, float).reshape(shape)); tot = tot + np.real(np.fft.ifftn(np.abs(Fh) ** 2))
                    tot = np.rint(tot).astype(np.int64).flatten()
                    first = (n, int(tot[0]), bool(np.all(tot[1:] == 0)), [x.tolist() for x in (X, X2, Y, Y2)])
        if stop_at and total >= stop_at: break
    return (total, first, ncls), bC

if __name__ == '__main__':
    maxC, log = int(sys.argv[1]), sys.argv[2]
    done = set()
    if os.path.exists(log):
        for line in open(log):
            if line.startswith('CASE '): done.add(line.split(' | ')[0])
    out = open(log, 'a')
    for spec in sys.argv[3:]:
        G = Group(spec); t0 = time.time(); nc = 0
        for H, gens in sorted(G.subgroups().items(), key=lambda kv: -len(kv[0])):
            key = 'CASE %s H=%s' % (spec, gens)
            if key in done: continue
            OC, _ = G.orbits([G.perm(h) for h in gens])
            if len(OC) > maxC or len(OC) < 3: continue
            t1 = time.time()
            res, bC = gs_dtied(G, gens, maxC)
            line = '%s | |H|=%d bC=%d res=%s t=%.0f' % (key, len(H), bC, None if res is None else (res[0], res[2]), time.time() - t1)
            if res and res[0]:
                line += ' FIRST ' + json.dumps({'n': res[1][0], 'paf0': res[1][1], 'ok': res[1][2], 'seqs': res[1][3]})
                print('FOUND', key, res[1][:3], flush=True)
            out.write(line + '\n'); out.flush(); nc += 1
        out.write('SPECDONE %s cases=%d t=%.0f\n' % (spec, nc, time.time() - t0)); out.flush()
        print('specdone', spec, nc, flush=True)

"""Tied Goethals-Seidel census: four +-1 sequences on G with sum PAF = 0 off 0 (H(4|G|) by the GS array),
A, B invariant under <H, n>, C under H, D = C o n^-1.  Usage: gs_tied.py maxC maxAB log spec ..."""
import sys, os, time, json, numpy as np
from gengroup import *

def gs_tied(G, Hgens, maxC, maxAB, stop_at=1):
    m = G.m
    sigs = [G.perm(h) for h in Hgens]
    OC, _ = G.orbits(sigs)
    bC = len(OC)
    if bC > maxC: return None, bC
    eng = Engine(G, sigs)
    XC = all_signs(bC); PC = eng.paf(XC, OC)
    K = len(eng.key_reps); R = keyvec(K)
    z = eng.labk[0]
    target = np.zeros(K, dtype=np.int64); target[z] = 4 * m
    kT = int(target @ R)
    H = G.closure(Hgens)
    seen = set(); total = 0; first = None; ncls = 0
    for n in G.mult_elems():
        cls = frozenset(G.mcomp(n, h) for h in H)
        if cls in seen or n in H: continue
        seen.add(cls); ncls += 1
        OA, _ = G.orbits(sigs + [G.perm(n)])
        bA = len(OA)
        if 2 * bA > maxAB: continue
        XA = all_signs(bA); PA = eng.paf(XA, OA)
        kA = lkeys(PA, R)
        kAB = (kA[:, None] + kA[None, :]).reshape(-1)
        ninv = tuple(f.minv(x) for f, x in zip(G.f, n))
        kp = eng.key_perm(G.perm(ninv))
        Rd = np.zeros_like(R); np.add.at(Rd, kp, R)
        kS = lkeys(PC, R + Rd)
        for i, j in match(kAB, kT - kS):
            a, b = divmod(i, len(XA))
            if np.array_equal(PA[a].astype(np.int64) + PA[b] + PC[j] + PC[j][kp], target):
                total += 1
                if first is None:
                    A = XA[a] @ OA; B = XA[b] @ OA; Cc = XC[j] @ OC
                    P = G.perm(n); D = np.zeros(m, dtype=np.int64); D[P] = Cc
                    tot = 0
                    shape = tuple(s for f in G.f for s in f.shape)
                    for f in (A, B, Cc, D):
                        Fh = np.fft.fftn(np.asarray(f, float).reshape(shape)); tot = tot + np.real(np.fft.ifftn(np.abs(Fh) ** 2))
                    tot = np.rint(tot).astype(np.int64).flatten()
                    first = (n, int(tot[0]), bool(np.all(tot[1:] == 0)), [x.tolist() for x in (A, B, Cc, D)])
        if stop_at and total >= stop_at: break
    return (total, first, ncls), bC

if __name__ == '__main__':
    maxC, maxAB, log = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    done = set()
    if os.path.exists(log):
        for line in open(log):
            if line.startswith('CASE '): done.add(line.split(' | ')[0])
    out = open(log, 'a')
    for spec in sys.argv[4:]:
        G = Group(spec); t0 = time.time(); nc = 0
        for H, gens in sorted(G.subgroups().items(), key=lambda kv: -len(kv[0])):
            key = 'CASE %s H=%s' % (spec, gens)
            if key in done: continue
            sigs = [G.perm(h) for h in gens]
            OC, _ = G.orbits(sigs)
            if len(OC) > maxC or len(OC) < 3: continue
            t1 = time.time()
            res, bC = gs_tied(G, gens, maxC, maxAB)
            line = '%s | |H|=%d bC=%d res=%s t=%.0f' % (key, len(H), bC, None if res is None else res[:1] + res[2:], time.time() - t1)
            if res and res[0]:
                line += ' FIRST ' + json.dumps({'n': res[1][0], 'paf0': res[1][1], 'ok': res[1][2], 'seqs': res[1][3]})
                print('FOUND', key, res[1][:3], flush=True)
            out.write(line + '\n'); out.flush(); nc += 1
        out.write('SPECDONE %s cases=%d t=%.0f\n' % (spec, nc, time.time() - t0)); out.flush()
        print('specdone', spec, nc, flush=True)

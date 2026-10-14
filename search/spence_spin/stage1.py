"""hier2 with (i) exact norm-consistent pair-count scoring of Z_13 pairs (DP over A groups and B columns), (ii) exact DP sampling of A candidates at a given energy
and of per-column energy allocations."""
import sys, time, subprocess, itertools
import numpy as np
from struct573 import build
from levels import make
import hashjoin

def main(seed, secs, npairs_scan, nrounds, nA_target, u=3, maxwall=170):
    S = build(u); LQ = make(S, 'Z2^2xZ13'); L13 = make(S, 'Z13'); G = S.G
    rng = np.random.default_rng(seed); N0 = LQ.norm; t_start = time.time()
    hits13 = []; s0 = seed
    while len(hits13) < npairs_scan:
        h, _, _, _ = hashjoin.run(S, L13, 100000, 100000, s0, chunk=100000); hits13 += h; s0 += 1000
    coordsQ = G.coords[LQ.reps]; tQ = coordsQ[:, 2] % 13
    t13 = G.coords[L13.reps][:, 2] % 13
    k13_of_t = {int(t): k for k, t in enumerate(t13)}
    colvars = [list(np.flatnonzero(tQ == t)) for t in range(13)]
    Cm = np.zeros((L13.RK, LQ.RK), dtype=np.int64)
    rep13 = [int(np.flatnonzero(L13.orbK == j)[0]) for j in range(L13.RK)]
    for j in range(L13.RK):
        for i in colvars[int(t13[rep13[j]])]: Cm[j, LQ.orbK[i]] += 1
    groups = [[R for R in range(LQ.RK) if Cm[j, R]] for j in range(L13.RK)]
    vs = [-11, -9, -1, 1, 9, 11]
    EMAX = N0 + 1
    def structs(a5, b13):
        gl = []
        for j in range(L13.RK):
            vars_ = groups[j]
            lst = [vals for vals in itertools.product(*[list(LQ.valsA[R].keys()) for R in vars_]) if sum(int(Cm[j, R]) * v for R, v in zip(vars_, vals)) == int(a5[j])]
            gl.append(np.array(lst, dtype=np.int64).reshape(len(lst), len(vars_)))
        bo = []
        for t in range(13):
            tg = int(b13[k13_of_t[t]]); bo.append(np.array([o for o in itertools.product(vs, repeat=4) if sum(o) == tg], dtype=np.int64).reshape(-1, 4))
        return gl, bo
    def ways(gl, bo):
        # A: energy contributions per group item: sum_R wA[R] * v^2 ; B: per column item energy
        wA = LQ.wA
        cA = np.zeros(EMAX); cA[0] = 1
        gen = []
        for j in range(L13.RK):
            e = ((gl[j] ** 2) * wA[groups[j]][None, :]).sum(1); gen.append(e)
            nc = np.zeros(EMAX)
            for ev in e:
                if ev < EMAX: nc[ev:] += cA[:EMAX - ev]
            cA = nc
        cB = np.zeros(EMAX); cB[0] = 1; cols = []
        for t in range(13):
            e = (bo[t] ** 2).sum(1); cols.append(e)
            nc = np.zeros(EMAX)
            for ev in e:
                if ev < EMAX: nc[ev:] += cB[:EMAX - ev]
            cB = nc
        return cA, cB, gen, cols
    scored = []
    for i, (a5, b13) in enumerate(hits13[:npairs_scan]):
        gl, bo = structs(a5, b13)
        if min(len(g) for g in gl) == 0 or min(len(o) for o in bo) == 0: continue
        cA, cB, _, _ = ways(gl, bo)
        tot = 0.0; bestE = None; bm = 0
        for eA in range(0, EMAX):
            if cA[eA] and (N0 - eA) % 3 == 0:
                eB = (N0 - eA) // 3
                if 0 <= eB < EMAX and cB[eB]:
                    w = cA[eA] * cB[eB]; tot += w
                    if w > bm: bm, bestE = w, eA
        if tot > 0: scored.append((tot, i, bestE))
    scored.sort(reverse=True)
    print('scanned %d Z13 pairs; with norm-consistent extensions: %d; best scores (log10 pairs): %s' % (npairs_scan, len(scored), [round(np.log10(s[0]), 1) for s in scored[:5]]), flush=True)
    best = None; tot_pairs = 0.0; tot_probes = 0.0
    for (tot, i, bestE) in scored[:3]:
        a5, b13 = hits13[i]
        gl, bo = structs(a5, b13)
        cA, cB, gen, cols = ways(gl, bo)
        # candidate energy splits ranked by cA*cB
        splits = sorted([(cA[eA] * cB[(N0 - eA) // 3], eA) for eA in range(EMAX) if cA[eA] and (N0 - eA) % 3 == 0 and cB[(N0 - eA) // 3]], reverse=True)
        for rnd in range(nrounds):
            if time.time() - t_start > maxwall: break
            w, EAs = splits[min(rnd, len(splits) - 1)]; EBs = (N0 - EAs) // 3
            # A sampling at exact energy EAs via group DP (suffix ways)
            suf = [None] * (L13.RK + 1); suf[L13.RK] = np.zeros(EMAX); suf[L13.RK][0] = 1
            for j in range(L13.RK - 1, -1, -1):
                nc = np.zeros(EMAX)
                for ev in gen[j]:
                    if ev < EMAX: nc[ev:] += suf[j + 1][:EMAX - ev]
                suf[j] = nc
            nA = nA_target
            Aarr = np.zeros((nA, LQ.RK), dtype=np.int16); rem = np.full(nA, EAs)
            for j in range(L13.RK):
                e = gen[j]
                # probability of item k: suf[j+1][rem - e_k]
                P = np.stack([np.where(rem - e[k] >= 0, suf[j + 1][np.clip(rem - e[k], 0, EMAX - 1)], 0.0) for k in range(len(e))], 1)
                P = P / P.sum(1, keepdims=True); cum = np.cumsum(P, 1); r = rng.random(nA)[:, None]
                k = (r > cum).sum(1); k = np.minimum(k, len(e) - 1)
                for col, R in enumerate(groups[j]): Aarr[:, R] = gl[j][k, col]
                rem = rem - e[k]
            assert np.all(rem == 0)
            Aarr = np.unique(Aarr, axis=0)
            HA = hashjoin.hashkey(hashjoin.keys_of(Aarr, LQ.MA))
            # per-column energy allocation sampled by DP
            sufB = [None] * 14; sufB[13] = np.zeros(EMAX); sufB[13][0] = 1
            for t in range(12, -1, -1):
                nc = np.zeros(EMAX)
                for ev in cols[t]:
                    if ev < EMAX: nc[ev:] += sufB[t + 1][:EMAX - ev]
                sufB[t] = nc
            remB = EBs; ecls = []
            for t in range(13):
                ev = np.unique(cols[t]); pw = np.array([(cols[t] == x).sum() * sufB[t + 1][remB - x] if remB - x >= 0 else 0 for x in ev], dtype=float)
                x = rng.choice(ev, p=pw / pw.sum()); ecls.append(int(x)); remB -= int(x)
            assert remB == 0
            ropts = [bo[t][(bo[t] ** 2).sum(1) == ecls[t]] for t in range(13)]
            fn = 'fill3_in.txt'
            with open(fn, 'w') as f:
                f.write('%d %d %d\n' % (LQ.RH, LQ.D, 13))
                for t in range(13):
                    f.write(' '.join(map(str, colvars[t])) + ' %d\n' % len(ropts[t]))
                    for o in ropts[t]: f.write(' '.join(map(str, o)) + '\n')
                for cc in range(LQ.D):
                    M = LQ.MB[cc] + LQ.MB[cc].T
                    for ii in range(LQ.RH): f.write(' '.join(map(str, M[ii])) + '\n')
                f.write('%d\n%d\n' % (LQ.target, len(Aarr)))
                f.write(' '.join(str(int(h) & 0xFFFFFFFFFFFFFFFF) for h in HA) + '\n')
            t0 = time.time()
            out = subprocess.run(['nice', '-n', '10', './fillB', fn, str(secs), str(seed + 31 * i + rnd), '30'], capture_output=True, text=True).stdout.strip().split('\n')
            pr = int(out[-1].split()[2]); tot_probes += pr; tot_pairs += pr * len(Aarr)
            dic = {}
            for ii, h in enumerate(HA): dic.setdefault(int(h) & 0xFFFFFFFFFFFFFFFF, []).append(ii)
            nver = 0
            for line in out:
                if not line.startswith('HIT'): continue
                parts = line.split(); h = int(parts[2]); bvals = np.array(list(map(int, parts[4:])), dtype=np.int64)
                b = np.zeros(LQ.RH, dtype=np.int64); pos = 0
                for t in range(13):
                    for k in range(4): b[colvars[t][k]] = bvals[pos]; pos += 1
                for ii in dic.get(h, []):
                    a = Aarr[ii].astype(np.int64)
                    if np.array_equal(hashjoin.keys_of(a[None], LQ.MA)[0] + hashjoin.keys_of(b[None], LQ.MB)[0], np.full(LQ.D, LQ.target)) and hashjoin.check_hit(LQ, a, b):
                        nver += 1
                        if True:
                            best = (a.copy(), b.copy()); _k = globals().setdefault('_nsaved', [0]); _k[0] += 1; np.save('wit_u%d_s%d_%d.npy' % (u, seed, _k[0]), np.concatenate([a, b]))
                            print('   *** VERIFIED Q-level solution at Z2^2xZ13:\n       a=%s\n       b=%s' % (a.tolist(), b.tolist()), flush=True)
            print('pair#%d (log10 norm-consistent pairs %.1f) split E_A=%d E_B=%d: |A table|=%d, B opts/col %s, probes %.3g, pairs %.3g, verified hits %d (%.0fs)' % (i, np.log10(tot), EAs, EBs, len(Aarr), [len(o) for o in ropts], pr, pr * len(Aarr), nver, time.time() - t0), flush=True)
            if best is not None: break
        if best is not None or time.time() - t_start > maxwall: break
    print('TOTAL pairs examined %.3g probes %.3g wall %.0fs' % (tot_pairs, tot_probes, time.time() - t_start))

if __name__ == '__main__':
    main(int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(float(sys.argv[5])))

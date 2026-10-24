"""spinD_drv.py -- d = 2 and d = 4 block-symmetric Spence structures over Z_v (v odd): H-invariant X (and Y) with blocks
   d = 4:  (X, mu X, mu^2 X, mu^3 X),  K = <H, mu>, K/H cyclic of order 4:   sum_{j<4} PAF_X(mu^j g) = -4  (g != 0)
   d = 2:  (X, mu X, Y, mu Y),         K/H cyclic of order 2:                 sum_{j<2} (PAF_X + PAF_Y)(mu^j g) = -4
all row sums +1 (gauge).  Exact: Gray-code enumeration of X (spinD.c), d = 2: self-join on the list.
usage: spinD_drv.py v d OHmax [OHmin]      selftest: spinD_drv.py selftest"""
import sys, os, math, subprocess, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'plain'))
import plain_spence3 as ps
import spin_drv as sd

def setupD(v, H, mu):
    K = ps.closure(list(H) + [mu], v)
    labH, repsH = ps.orbits(v, H); labK, repsK = ps.orbits(v, K)
    OH = len(repsH); OK = len(repsK)
    sizeH = np.bincount(labH, minlength=OH)
    negl = labK[(-np.array(repsK)) % v]
    seen = set(); cls = []
    for a in range(1, OK):
        if a in seen: continue
        seen.add(a); seen.add(int(negl[a])); cls.append(a)
    return labH, OH, sizeH, [repsK[a] for a in cls], OK

def build_inputD(v, H, mu, d, path):
    labH, OH, sizeH, gs, OK = setupD(v, H, mu)
    C = len(gs)
    SB = np.zeros((C, OH, OH), dtype=np.int64)
    for c, g in enumerate(gs):
        for j in range(d):
            SB[c] += sd.pair_forms(v, labH, OH, (pow(mu, j, v) * g) % v)
    with open(path, 'w') as f:
        f.write(f"{1 if d == 4 else 2} {C} {OH}\n")
        f.write(' '.join(map(str, sizeH.tolist())) + '\n')
        f.write(' '.join(map(str, SB.ravel().tolist())) + '\n')
        f.write(' '.join(['-4'] * C) + '\n')
    return dict(C=C, OH=OH, OK=OK, labH=labH, sizeH=sizeH)

def seqs_of(v, info, d, mu, gray, gray2=None):
    OH = info['OH']
    def xs(gr): return np.array([1 - 2 * ((gr >> i) & 1) for i in range(OH)], dtype=np.int64)[info['labH']]
    X = xs(gray); h = np.arange(v)
    if d == 4: return [X[(pow(mu, j, v) * h) % v] for j in range(4)]
    Y = xs(gray2); return [X, X[(mu * h) % v], Y, Y[(mu * h) % v]]

def verify(v, seqs):
    tot = np.zeros(v, dtype=np.int64); h = np.arange(v)
    for s in seqs:
        for g in range(v): tot[g] += int(s @ s[(h + g) % v])
    return bool(tot[0] == 4 * v and np.all(tot[1:] == -4) and all(int(s.sum()) == 1 for s in seqs))

def run_pairD(v, H, mu, d, tag=None):
    path = os.path.join(HERE, f'tmpD_{os.getpid()}.txt')
    labH, OH, sizeH, gs, OK = setupD(v, H, mu)
    if not sd.rho_feasible(sizeH): return 0, [], dict(OH=OH, OK=OK, C=len(gs), L=0)
    info = build_inputD(v, H, mu, d, path)
    r = subprocess.run([os.path.join(HERE, 'spinD'), path], capture_output=True, text=True)
    count = 0; sols = []; L = None
    for line in r.stdout.splitlines():
        if line.startswith('COUNT'): count = int(line.split()[1])
        elif line.startswith('LIST'): L = int(line.split()[1])
        elif line.startswith('SOL'): sols.append(tuple(int(t) for t in line.split()[1:]))
    info['L'] = L; info['raw'] = r.stdout[-200:]
    return count, sols, info

def brute(v, H, mu, d):
    labH, OH, sizeH, gs, OK = setupD(v, H, mu)
    X = ps.enum_list(OH, sizeH, rowsum=1)
    h = np.arange(v)
    W = []
    for x in X:
        s = x.astype(np.int64)[labH]
        w = []
        for g in gs:
            t = 0
            for j in range(d):
                gg = (pow(mu, j, v) * g) % v
                t += int(s @ s[(h + gg) % v])
            w.append(t)
        W.append(w)
    W = np.array(W, dtype=np.int64).reshape(len(X), len(gs))
    if d == 4: return int(np.sum(np.all(W == -4, axis=1)))
    cnt = 0; tgt = -4
    from collections import Counter
    cc = Counter(map(tuple, W.tolist()))
    for w, m in cc.items():
        o = tuple(tgt - a for a in w)
        cnt += m * cc.get(o, 0)
    return cnt

def pairsD(v, d, OHmax, OHmin=2):
    subs = ps.all_subgroups(v); U = ps.units(v); out = {}
    for H in subs:
        labH, repsH = ps.orbits(v, H); OH = len(repsH)
        if not (OHmin <= OH <= OHmax): continue
        Hs = set(H); gH = sd.small_gens(H, v)
        for mu in U:
            if mu in Hs: continue
            if pow(mu, d, v) not in Hs: continue
            if any(pow(mu, j, v) in Hs for j in range(1, d)): continue
            K = ps.closure(gH + [mu], v)
            out.setdefault((H, K), mu)
    return [(H, mu) for (H, K), mu in out.items()]

def selftest():
    bad = 0; nz = 0; n = 0
    for v in (9, 13, 15, 17, 21, 25, 27, 29, 31, 33, 35, 37, 39, 41, 45, 49, 51, 55, 57, 63):
        for d in (2, 4):
            for H, mu in pairsD(v, d, 15, 2):
                labH, OH, sizeH, gs, OK = setupD(v, H, mu)
                if OH > 13: continue
                c, sols, info = run_pairD(v, H, mu, d)
                b = brute(v, H, mu, d)
                n += 1
                if c != b: bad += 1; print('MISMATCH', v, d, len(H), mu, c, b)
                if b: nz += 1; print('nonzero', v, d, 'OH', OH, 'h', len(H), 'mu', mu, 'count', c, 'brute', b, flush=True)
                if c and d == 4:
                    for s in sols[:2]:
                        ok = verify(v, seqs_of(v, info, d, mu, s[0])); print('   verify', ok)
                if c and d == 2:
                    okc = 0
                    for s in sols[:4]:
                        okc += verify(v, seqs_of(v, info, d, mu, s[0], s[1]))
                    print('   verify', okc, 'of', len(sols[:4]))
    print('selftest: pairs', n, 'mismatches', bad, 'with nonzero count', nz)

if __name__ == '__main__':
    if sys.argv[1] == 'selftest': selftest(); sys.exit()
    v = int(sys.argv[1]); d = int(sys.argv[2]); OHmax = int(sys.argv[3]); OHmin = int(sys.argv[4]) if len(sys.argv) > 4 else 2
    prs = pairsD(v, d, OHmax, OHmin)
    print(f'v = {v}, d = {d}: {len(prs)} pairs (H, mu) with {OHmin} <= O_H <= {OHmax}', flush=True)
    for H, mu in sorted(prs, key=lambda hm: (len(ps.orbits(v, hm[0])[1]), -len(hm[0]))):
        t0 = time.time()
        try: count, sols, info = run_pairD(v, H, mu, d)
        except Exception as ex: print('ERR', len(H), mu, repr(ex)[:100], flush=True); continue
        rec = dict(v=v, d=d, OH=info['OH'], h=len(H), mu=mu, OK=info['OK'], C=info['C'], L=info.get('L'), count=count, secs=round(time.time() - t0, 1))
        if count:
            okv = []
            for s in sols[:20]:
                sq = seqs_of(v, info, d, mu, s[0], s[1] if len(s) > 1 else None)
                if verify(v, sq): okv.append([q.tolist() for q in sq]); break
            rec['verified'] = bool(okv); rec['sol'] = okv[0] if okv else None
            print('FOUND', rec, flush=True)
        print({k: val for k, val in rec.items() if k != 'sol'}, flush=True)

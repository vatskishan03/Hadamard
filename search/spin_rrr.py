"""RRR for spin structures, vectorised, with restarts and seeds.  Usage: spin_rrr.py v hgen mu iters_per_restart restarts beta seed"""
import sys, math, time
import numpy as np

def orbits(v, H):
    seen = np.full(v, -1); reps = []
    for x in range(v):
        if seen[x] >= 0: continue
        o = sorted({h * x % v for h in H})
        for y in o: seen[y] = len(reps)
        reps.append(o)
    return reps, seen

def closure(gens, v):
    S = {1}
    while True:
        T = S | {s * g % v for s in S for g in gens}
        if T == S: return sorted(S)
        S = T

def setup(v, hg, mu):
    H = closure([hg], v); K = closure([hg, mu], v)
    assert len(K) == 3 * len(H) and (v - 1) not in K
    OH, segH = orbits(v, H); OK_, segK = orbits(v, K)
    nH, nK = len(OH), len(OK_)
    b = int(round(math.sqrt((4 * v - 1) / 3))); assert 1 + 3 * b * b == 4 * v
    hs = len(H)
    s_plus = (b - 1) // hs if (b - 1) % hs == 0 else (-b - 1) // hs
    assert (1 + hs * s_plus) in (b, -b)
    m = nH - 1; npos = (m + s_plus) // 2
    w = np.exp(2j * np.pi * np.outer(np.arange(v), np.arange(v)) / v)
    reps = [o[0] for o in OH]
    eta = np.zeros((nH, nH), dtype=complex)
    for i, o in enumerate(OH):
        eta[:, i] = w[np.ix_(reps, o)].sum(axis=1)
    L = eta[1:, 1:]
    jmu = np.array([segH[reps[j] * mu % v] for j in range(nH)])
    seenj = np.zeros(nH, bool); triples = []
    for j in range(1, nH):
        if seenj[j]: continue
        t = [j, jmu[j], jmu[jmu[j]]]
        for x in t: seenj[x] = True
        triples.append(t)
    tri = np.array(triples) - 1
    repK = [o[0] for o in OK_]
    negK = np.array([segK[(-repK[q]) % v] for q in range(nK)])
    pairs = []; seenq = np.zeros(nK, bool); seenq[0] = True
    for q in range(1, nK):
        if seenq[q]: continue
        pairs.append((q, negK[q])); seenq[q] = seenq[negK[q]] = True
    NP = len(pairs)
    zeta = np.zeros((nK, nK), dtype=complex)
    for q, o in enumerate(OK_):
        zeta[:, q] = w[np.ix_(repK, o)].sum(axis=1)
    tri_kq = [segK[reps[t[0]]] for t in triples]
    allA = np.array(np.meshgrid(*[[1, -1]] * NP)).reshape(NP, -1).T
    avec = np.zeros((len(allA), nK)); avec[:, 0] = 1
    for p, (q1, q2) in enumerate(pairs):
        avec[:, q1] = allA[:, p]; avec[:, q2] = -allA[:, p]
    A2 = np.abs((avec @ zeta.T)[:, tri_kq]) ** 2
    T = 4 * v - A2
    ok = np.all(T >= -1e-9, axis=1)
    return dict(v=v, m=m, npos=npos, L=L, tri=tri, T=np.maximum(T[ok], 0), allA=allA[ok], hs=hs, nH=nH, NP=NP)

def run(S, iters, restarts, beta, seed, verbose=True):
    rng = np.random.default_rng(seed)
    m, npos, L, tri, T, allA = S['m'], S['npos'], S['L'], S['tri'], S['T'], S['allA']
    sqT = np.sqrt(T)
    Lr = np.vstack([L.real, L.imag]); Linv = np.linalg.pinv(Lr)
    def P1(x):
        y = np.full(m, -1.0); y[np.argpartition(-x, npos)[:npos]] = 1.0; return y
    def P2(x):
        be = 1 + L @ x
        mag2 = np.abs(be) ** 2
        nrm = np.sqrt(mag2[tri].sum(axis=1))
        d = ((sqT - nrm[None, :]) ** 2).sum(axis=1)
        a = int(np.argmin(d))
        sc = sqT[a] / np.maximum(nrm, 1e-9)
        fac = np.empty(m); fac[tri.ravel()] = np.repeat(sc, 3)
        nb = be * fac - 1
        return Linv @ np.concatenate([nb.real, nb.imag])
    def resid(y):
        be = 1 + L @ y
        s = (np.abs(be) ** 2)[tri].sum(axis=1)
        d = ((T - s[None, :]) ** 2).sum(axis=1)
        a = int(np.argmin(d)); return d[a], a
    t0 = time.time()
    for r in range(restarts):
        x = rng.standard_normal(m) * 3
        best = None
        for it in range(iters):
            p1 = P1(x)
            x = x + beta * (P2(2 * p1 - x) - p1)
            if it % 25 == 0:
                rr, a = resid(p1)
                if best is None or rr < best: best = rr
                if rr < 1e-6:
                    print(f'SOLUTION restart {r} iter {it} [{time.time()-t0:.0f}s]: A {allA[a].tolist()} B {[1] + p1.astype(int).tolist()}', flush=True)
                    return True
        if verbose:
            print(f'restart {r}: best {best:.0f} [{time.time()-t0:.0f}s]', flush=True)
    print(f'no solution in {restarts} restarts [{time.time()-t0:.0f}s]', flush=True)
    return False

if __name__ == '__main__':
    v, hg, mu = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    iters, restarts, beta, seed = int(sys.argv[4]), int(sys.argv[5]), float(sys.argv[6]), int(sys.argv[7])
    S = setup(v, hg, mu)
    print(f'v={v} |H|={S["hs"]} orbits {S["nH"]} triples {len(S["tri"])} A pairs {S["NP"]} admissible A {len(S["T"])} npos {S["npos"]}/{S["m"]}', flush=True)
    run(S, iters, restarts, beta, seed)

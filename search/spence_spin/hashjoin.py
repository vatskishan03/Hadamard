"""Random 2-list hash join at a quotient level: sample compressed A (row sum 2) and B (row sum 0) from the exact sign-level pools, key = class-PAF vectors,
look for T^A + T^B = target.  Hits are certificates of non-emptiness of the level (verified independently by direct PAF on Q)."""
import numpy as np, sys, time, math
from struct573 import build
from levels import make

PR = np.array([0x9E3779B97F4A7C15, 0xC2B2AE3D27D4EB4F, 0x165667B19E3779F9, 0xD6E8FEB86659FD93, 0xFF51AFD7ED558CCD, 0xC4CEB9FE1A85EC53,
               0x2545F4914F6CDD1D, 0x94D049BB133111EB, 0xBF58476D1CE4E5B9, 0x9FB21C651E98DF25, 0x1B873593CC9E2D51, 0xE6546B64AB4F5C0F], dtype=np.uint64).view(np.int64)

def compress_matrices(S, L):
    cntA = np.zeros((S.kK, L.nQ)); np.add.at(cntA, (S.orbK, L.cos), 1)
    cntB = np.zeros((S.kH, L.nQ)); np.add.at(cntB, (S.orbH, L.cos), 1)
    repA = [int(np.flatnonzero(L.orbK == r)[0]) for r in range(L.RK)]; repB = [int(np.flatnonzero(L.orbH == r)[0]) for r in range(L.RH)]
    return cntA[:, repA], cntB[:, repB]

def keys_of(V, M):
    V = V.astype(np.float64)
    return np.stack([((V @ M[c].astype(np.float64)) * V).sum(1) for c in range(M.shape[0])], 1).round().astype(np.int64)

def hashkey(K):
    with np.errstate(over='ignore'):
        h = np.zeros(len(K), dtype=np.int64)
        for c in range(K.shape[1]): h = h + K[:, c] * PR[c % len(PR)]
        return h

def run(S, L, nA, nB, seed, chunk=500000):
    rng = np.random.default_rng(seed)
    CA, CB = compress_matrices(S, L)
    t0 = time.time()
    hs, Vs = [], []
    for s in range(0, nA, chunk):
        m = min(chunk, nA - s)
        XA = S.sample('A', 2, m, rng); VA = np.rint(XA @ CA).astype(np.int16)
        KA = keys_of(VA, L.MA); hs.append(hashkey(KA)); Vs.append(VA)
    HA = np.concatenate(hs); VA = np.concatenate(Vs)
    order = np.argsort(HA); HA = HA[order]; VA = VA[order]
    t1 = time.time()
    hits = []; hits_total = [0]
    for s in range(0, nB, chunk):
        m = min(chunk, nB - s)
        XB = S.sample('B', 0, m, rng); VB = np.rint(XB @ CB).astype(np.int16)
        KB = keys_of(VB, L.MB)
        need = hashkey(L.target - KB)
        pos = np.searchsorted(HA, need); pos[pos >= len(HA)] = len(HA) - 1
        ok = HA[pos] == need
        nh = int(ok.sum()); hits_total[0] += nh
        for j in np.flatnonzero(ok)[:50]:
            # collision check by exact key
            i = pos[j]
            while i < len(HA) and HA[i] == need[j]:
                a = VA[i].astype(np.int64); b = VB[j].astype(np.int64)
                if np.array_equal(keys_of(a[None], L.MA)[0] + keys_of(b[None], L.MB)[0], np.full(L.D, L.target)): hits.append((a, b))
                i += 1
    return hits, t1 - t0, time.time() - t1, hits_total[0]

def check_hit(L, a, b):
    abar = a[L.orbK]; bbar = b[L.orbH]
    def paf(v, q): return int(v @ v[L.Qadd[:, q]])
    ok = all(paf(abar, q) + paf(bbar, q) + paf(bbar[L.pmu], q) + paf(bbar[L.pmu2], q) == L.target for q in range(1, L.nQ))
    return ok and int(abar.sum()) == 2 and int(bbar.sum()) == 0

if __name__ == '__main__':
    u, name, nA, nB, seed = int(sys.argv[1]), sys.argv[2], int(float(sys.argv[3])), int(float(sys.argv[4])), int(sys.argv[5])
    S = build(u); L = make(S, name)
    hits, ta, tb, ht = run(S, L, nA, nB, seed)
    print(name, 'u=%d' % u, 'nA=%d nB=%d raw hash hits=%d  (build %.1fs, probe %.1fs)  density %.3g' % (nA, nB, ht, ta, tb, ht / (nA * nB)), flush=True)
    for a, b in hits[:3]:
        print('  verified:', check_hit(L, a, b), ' a=', a.tolist(), ' b=', b.tolist())
    if hits:
        np.save('wit_%s_u%d_hj.npy' % (name.replace('^', ''), u), np.concatenate([hits[0][0], hits[0][1]]))

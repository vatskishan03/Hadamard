"""Rebuild a spin solution (A on K-orbit pairs, B on H-orbits) as sequences over Z_v and check the GS skew Hadamard."""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'verify'))
import verify
def closure(gens, v):
    S = {1}
    while True:
        T = S | {s * g % v for s in S for g in gens}
        if T == S: return sorted(S)
        S = T
def orbits(v, H):
    seen = np.full(v, -1); reps = []
    for x in range(v):
        if seen[x] >= 0: continue
        o = sorted({h * x % v for h in H})
        for y in o: seen[y] = len(reps)
        reps.append(o)
    return reps, seen
def build(v, hg, mu, Avals, Bvals):
    H = closure([hg], v); K = closure([hg, mu], v)
    OH, segH = orbits(v, H); OK_, segK = orbits(v, K)
    repK = [o[0] for o in OK_]
    negK = [segK[(-repK[q]) % v] for q in range(len(OK_))]
    pairs = []; seen = [False] * len(OK_); seen[0] = True
    for q in range(1, len(OK_)):
        if seen[q]: continue
        pairs.append((q, negK[q])); seen[q] = seen[negK[q]] = True
    aK = np.zeros(len(OK_), dtype=np.int64); aK[0] = 1
    for p, (q1, q2) in enumerate(pairs):
        aK[q1] = Avals[p]; aK[q2] = -Avals[p]
    a = aK[segK]
    bH = np.array(Bvals, dtype=np.int64)
    b = bH[segH]
    t = np.arange(v)
    c = b[t * mu % v]; d = b[t * mu * mu % v]
    return [a, b, c, d]
if __name__ == '__main__':
    v, hg, mu = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    A = json.loads(sys.argv[4]); B = json.loads(sys.argv[5])
    seqs = build(v, hg, mu, A, B)
    tot = verify.paf_sum(seqs)
    print('PAF sum at 0:', tot[0], 'off-peak nonzero:', int(np.count_nonzero(tot[1:])))
    Hm = verify.goethals_seidel(*map(verify.circulant, seqs)).astype(np.int64)
    n = 4 * v
    print('Hadamard', np.array_equal(Hm @ Hm.T, n * np.eye(n, dtype=np.int64)), 'skew', np.array_equal(Hm + Hm.T, 2 * np.eye(n, dtype=np.int64)))

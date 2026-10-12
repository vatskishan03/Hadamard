"""Exhaustive skew Goethals-Seidel search on a multiplier structure: A skew-type (A(0)=1, A(-x)=-A(x)),
B, C, D constant on H-orbits, all admissible row-sum patterns, 4-list exact matcher find4b.
Each solution is verified: GS array H with H H^T = 4n I and H + H^T = 2 I (skew-Hadamard)."""
import os, struct, subprocess, sys, tempfile, time, json, math
import numpy as np
sys.path.insert(0, '/home/user/Hadamard/attack'); sys.path.insert(0, '/home/user/Hadamard/verify')
import spence_lib as L, count_exact as X
import verify
from skew_ktree import pool_exact
HERE = '/home/user/Hadamard/attack'

def run(n, gens, cap=4e10, K=2, verbose=True, only=None):
    G = L.AbGroup([n]); st = L.Structure(G, [G.perm_from_matrix(np.array([[g]])) for g in gens])
    rep = np.array([np.flatnonzero(st.orb == o)[0] for o in range(st.O)])
    partner = st.orb[G.neg()[rep]]
    assert not np.any(partner[1:] == np.arange(1, st.O)), "-1 fixes an orbit"
    pairs = [(o, partner[o]) for o in range(1, st.O) if o < partner[o]]
    idx = np.arange(2 ** len(pairs))
    s = 1 - 2 * ((idx[:, None] >> np.arange(len(pairs))) & 1)
    YA = np.zeros((len(idx), st.O), dtype=np.int64); YA[:, 0] = 1
    for j, (o, q) in enumerate(pairs): YA[:, o] = s[:, j]; YA[:, q] = -s[:, j]
    pm = n % 4
    VA = np.rint((st.class_values(YA.astype(float)) - pm) / 4).astype(np.int16)
    lim = int(math.isqrt(4 * n)); pats = []
    for b in range(1, lim + 1, 2):
        for c in range(b, lim + 1, 2):
            d2 = 4 * n - 1 - b * b - c * c
            if d2 < c * c: continue
            d = math.isqrt(d2)
            if d * d == d2 and all(st.pool(r) > 0 for r in (b, c, d)): pats.append((b, c, d))
    if only: pats = [p for p in pats if p in only]
    if verbose: print(f"n={n} gens={gens} orbits={st.O} classes={st.C} skew pool={len(YA)} patterns={pats}", flush=True)
    cache = {}
    def vecs(r):
        if r not in cache:
            Y = pool_exact(st, r).astype(np.float64)
            V = np.empty((len(Y), st.C), np.int16)
            step = max(1000, int(1.5e8 // (st.C * st.O)))
            for s0 in range(0, len(Y), step):
                V[s0:s0 + step] = np.rint((st.class_values(Y[s0:s0 + step]) - pm) / 4)
            cache[r] = (Y.astype(np.int8), V)
        return cache[r]
    found = []
    for q in pats:
        sizes = [len(YA)] + [int(st.pool(r)) for r in q]
        work = sizes[0] * sizes[3] + sizes[1] * sizes[2]
        if work > cap:
            if verbose: print(f"  pattern {q}: pools {sizes} work {work:.2e} > cap", flush=True)
            continue
        Ls = [(YA.astype(np.int8), VA)] + [vecs(r) for r in q]
        fd, path = tempfile.mkstemp(suffix='.bin', dir=HERE); os.close(fd)
        try:
            with open(path, 'wb') as f:
                f.write(struct.pack('iiiii', *sizes, st.C)); f.write(np.full(st.C, -pm, dtype=np.int16).tobytes())
                for l in Ls: f.write(np.ascontiguousarray(l[1]).tobytes())
            t = time.time()
            out = subprocess.run([os.path.join(HERE, 'find4b'), path, str(K)], capture_output=True, text=True, check=True).stdout
        finally:
            os.remove(path)
        total = int(out.strip().split('TOTAL')[-1])
        if verbose: print(f"  pattern {q}: pools {sizes} work {work:.2e}: {total} solutions [{time.time()-t:.1f}s]", flush=True)
        for line in out.splitlines():
            if not line.startswith('SOL'): continue
            ids = list(map(int, line.split()[1:]))
            seqs = [Ls[k][0][ids[k]].astype(np.int64)[st.orb] for k in range(4)]
            tot = verify.paf_sum(seqs)
            if not (tot[0] == 4 * n and not tot[1:].any()): continue
            H = verify.goethals_seidel(*map(verify.circulant, seqs)).astype(np.float64)   # BLAS; exact below 2^53
            had = np.array_equal(H @ H.T, 4 * n * np.eye(4 * n))
            skew = np.array_equal(H + H.T, 2 * np.eye(4 * n))
            found.append(dict(n=n, gens=list(gens), pattern=list(q), hadamard=bool(had), skew=bool(skew),
                              seqs=[''.join('+' if x > 0 else '-' for x in sq) for sq in seqs]))
            if verbose: print(f"    VERIFIED: Hadamard={had} skew={skew}", flush=True)
    return found

if __name__ == '__main__':
    n = int(sys.argv[1]); gens = [int(g) for g in sys.argv[2].split(',')]
    cap = float(sys.argv[3]) if len(sys.argv) > 3 else 4e10
    only = [tuple(int(x) for x in p.split(',')) for p in sys.argv[4:]] or None
    found = run(n, gens, cap=cap, only=only)
    if found:
        with open(f'{HERE}/skew_hits.json', 'a') as f:
            for h in found: f.write(json.dumps(h) + '\n')
    print(f"done: {len(found)} verified")

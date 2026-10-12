"""Exact sweep of Goethals-Seidel multiplier structures for the unresolved Hadamard orders 4n.

The repository's earlier sweep (FRONTIER.md) searched every structure with at most 18 cosets and
left those with more than 22 orbits unsearched. The skew campaign (NOTES section 10) found
solutions in several structures that the Gaussian model put at 0.1-0.2. So this sweep ignores the
model.
- It takes every subgroup H of Z_n^* whose structure has MINO..MAXO orbits.
- It takes every row-sum pattern a^2 + b^2 + c^2 + d^2 = 4n with pair work N_1 N_4 + N_2 N_3 <= CAP.
- It drops patterns that fail the compression test at any prime d | n.
- It runs the exact four-list matcher find4b on the rest.

  python3 had_sweep.py list MINO MAXO CAP n ...     # candidates with work and compression verdict
  python3 had_sweep.py run  n gens a,b,c,d [CAP]    # exact search of one structure and pattern
"""
import itertools, json, math, os, struct, subprocess, sys, tempfile, time
from collections import defaultdict
import numpy as np
sys.path.insert(0, '/home/user/Hadamard/attack'); sys.path.insert(0, '/home/user/Hadamard/verify')
import spence_lib as L
from gs_sweep_subgroups import subgroups
import skew_compress as SC
from skew_ktree import pool_exact
import verify
HERE = '/home/user/Hadamard/attack'


def structure(n, gens):
    G = L.AbGroup([n])
    return L.Structure(G, [G.perm_from_matrix(np.array([[g]])) for g in gens])


def patterns(st, n):
    lim = int(math.isqrt(4 * n)); out = []
    for a in range(1, lim + 1, 2):
        for b in range(a, lim + 1, 2):
            for c in range(b, lim + 1, 2):
                d2 = 4 * n - a * a - b * b - c * c
                if d2 < c * c: continue
                d = math.isqrt(d2)
                if d * d == d2 and all(st.pool(r) > 0 for r in (a, b, c, d)): out.append((a, b, c, d))
    return out


def compression_ok(st, n, pat, d, cap=150_000):
    """False if the compressions to Z_d cannot satisfy the norm equations; None if too large to decide."""
    opts = [SC.compressed_options(st, n, d, r=r, cap=cap) for r in pat]
    if any(o is None for o in opts): return None
    Ns = [SC.norms(o, d) for o in opts]
    left = set()
    for a in Ns[0]:
        for b in Ns[1]: left.add(tuple(x + y for x, y in zip(a, b)))
    for c in Ns[2]:
        for e in Ns[3]:
            if tuple(4 * n - x - y for x, y in zip(c, e)) in left: return True
    return False


def primes(n):
    return [p for p in range(2, n) if n % p == 0 and all(p % q for q in range(2, int(p ** 0.5) + 1))]


def run(n, gens, pat, cap=1e11, K=2):
    st = structure(n, gens); pm = n % 4
    sizes = [int(st.pool(r)) for r in pat]
    order = np.argsort(sizes)                 # pair the smallest list with the largest
    roles = [pat[order[0]], pat[order[1]], pat[order[2]], pat[order[3]]]
    work = sizes[order[0]] * sizes[order[3]] + sizes[order[1]] * sizes[order[2]]
    print(f"n={n} gens={gens} orbits={st.O} classes={st.C} rows={pat} pools={sizes} work={work:.2e}", flush=True)
    if work > cap: return []
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
    # find4b pairs (list0, list3) on the left and (list1, list2) on the right
    Ls = [vecs(roles[0]), vecs(roles[1]), vecs(roles[2]), vecs(roles[3])]
    fd, path = tempfile.mkstemp(suffix='.bin', dir=HERE); os.close(fd)
    t = time.time()
    try:
        with open(path, 'wb') as f:
            f.write(struct.pack('iiiii', *[len(l[1]) for l in Ls], st.C)); f.write(np.full(st.C, -pm, dtype=np.int16).tobytes())
            for l in Ls: f.write(np.ascontiguousarray(l[1]).tobytes())
        out = subprocess.run([os.path.join(HERE, 'find4b'), path, str(K)], capture_output=True, text=True, check=True).stdout
    finally:
        os.remove(path)
    total = int(out.strip().split('TOTAL')[-1])
    print(f"  {total} solutions [{time.time() - t:.0f}s]", flush=True)
    found = []
    for line in out.splitlines():
        if not line.startswith('SOL'): continue
        ids = list(map(int, line.split()[1:]))
        seqs = [Ls[k][0][ids[k]].astype(np.int64)[st.orb] for k in range(4)]
        tot = verify.paf_sum(seqs)
        if not (tot[0] == 4 * n and not tot[1:].any()): continue
        H = verify.goethals_seidel(*map(verify.circulant, seqs)).astype(np.float64)
        had = np.array_equal(H @ H.T, 4 * n * np.eye(4 * n))
        print(f"  VERIFIED Hadamard={had}", flush=True)
        if had:
            rec = dict(n=n, gens=list(gens), pattern=[int(x.sum()) for x in seqs], hadamard=True,
                       seqs=[''.join('+' if x > 0 else '-' for x in s) for s in seqs])
            found.append(rec)
            with open(os.path.join(HERE, 'had_hits.json'), 'a') as f: f.write(json.dumps(rec) + '\n')
    return found


if __name__ == '__main__':
    if sys.argv[1] == 'list':
        mino, maxo, cap = int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
        for n in map(int, sys.argv[5:]):
            ds = [d for d in primes(n) if d < n]
            for H, gens in subgroups(n).items():
                O = len({min(t * h % n for h in H) for t in range(n)})
                if not (mino <= O <= maxo): continue
                st = structure(n, list(gens))
                for pat in patterns(st, n):
                    sizes = sorted(st.pool(r) for r in pat)
                    work = sizes[0] * sizes[3] + sizes[1] * sizes[2]
                    if work > cap: continue
                    verdict = 'ok'
                    for d in ds:
                        c = compression_ok(st, n, pat, d)
                        if c is False: verdict = f'DEAD(d={d})'; break
                        if c is None: verdict = 'ok?'
                    print(f"{work:.2e} n={n} 4n={4*n} gens={','.join(map(str, gens))} |H|={len(H)} O={st.O} C={st.C} "
                          f"rows={pat} {verdict}", flush=True)
    else:
        n = int(sys.argv[2]); gens = [int(x) for x in sys.argv[3].split(',')]
        pat = tuple(int(x) for x in sys.argv[4].split(','))
        cap = float(sys.argv[5]) if len(sys.argv) > 5 else 1e11
        run(n, gens, pat, cap)

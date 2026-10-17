"""Two-block searches over multiplier structures: D-optimal designs and Legendre pairs, on cyclic and non-cyclic groups.

Two +-1 functions A, B on an abelian group G of odd order v with
    PAF_A(s) + PAF_B(s) = T   for every s != 0
give
  * T = +2: a D-optimal design of order 2v, [[A, B], [-B^T, A^T]] with A A^T + B B^T = (2v - 2) I + 2 J
    (row sums with r_A^2 + r_B^2 = 4v - 2; group-developed blocks over any abelian group commute);
  * T = -2: a Legendre pair (row sums 1, 1), which gives a Hadamard matrix of order 2v + 2 (only cyclic G counts
    for the Legendre-pair problem itself);
  * T = 0, v even (`pg`): a periodic Golay pair (row sums with r_A^2 + r_B^2 = 2v), which gives Miyamoto-block
    Hadamard matrices (hadamard/MIYAMOTO_BLOCKS.md#periodic-golay-pairs) when v = (q - 1)/2.
Both sequences are constant on the orbits of a group H of automorphisms (multipliers of Z_v, or for v with a split
prime p^2 || v the GL(2, p) subgroups of propus_group.candidates_all times subgroups of Z_m^*). The autocorrelations
are constant on the classes of <H, -1>. For every structure with at most MAXO orbits and every row pattern the two
complete pools are listed (meet in the middle) and joined exactly on all classes (sorted 64-bit keys, then an exact
check). Larger pools (up to --psdmax) first go through the power-spectrum test |X^(chi)|^2 <= 2v - T of psd_pool.c,
which every member of a solution passes. Hits are verified as matrices and recorded in $HADAMARD_WORK/two_circ_hits.jsonl.

  python3 search/two_circ.py dopt v1,v2,.. [--maxo 34] [--maxpool 4e7] [--psdmax 3e10] [--group]
  python3 search/two_circ.py lp l1,l2,.. [--maxo 34] [--maxpool 4e7]
  python3 search/two_circ.py pg v1,v2,.. [--maxo 34] [--maxpool 4e7] [--psdmax 3e10]
"""
import argparse, json, math, sys, time
import numpy as np
from pathlib import Path as _P
sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
import spence_lib as L
import propus as PR
import propus_group as PG
from skew_ktree import pool_exact


def cyclic_structures(v, maxo):
    G = L.AbGroup([v])
    for S, gs in PR.subgroup_gens(v).items():
        st = L.Structure(G, [G.perm_from_matrix(np.array([[g % v]])) for g in gs])
        if st.O <= maxo: yield f"Z{v}<{','.join(map(str, gs))}>", G, [np.array([[g % v]]) for g in gs], st


def group_structures(v, maxo):
    for name, G, mats, st_ad, st_b in PG.structures(v, 1, maxo):
        if st_b.O <= maxo: yield name, G, mats, st_b


def patterns(v, T):
    if T == -2: return [(1, 1)]
    if T == 0:
        out = []
        for a in range(0, math.isqrt(2 * v) + 1, 2):
            b2 = 2 * v - a * a
            b = math.isqrt(b2)
            if b * b == b2 and b % 2 == 0 and a <= b: out.append((a, b))
        return out
    out = []; lim = math.isqrt(4 * v - 2)
    for a in range(1, lim + 1, 2):
        b2 = 4 * v - 2 - a * a
        b = math.isqrt(b2) if b2 > 0 else -1
        if b > 0 and b * b == b2 and b % 2 and a <= b: out.append((a, b))
    return out


def keys(V, w):
    return (V.astype(np.int64) * w[None, :]).sum(1)          # 64-bit wrap-around is fine for a hash


def join2(VA, VB, target):
    """All (i, j) with VA[i] + VB[j] = target (exact): sorted hash keys, then an exact check of the candidates."""
    rng = np.random.default_rng(12345)
    w = rng.integers(1, 2 ** 62, size=VA.shape[1], dtype=np.int64)
    ka = keys(target[None, :] - VA, w); kb = keys(VB, w)
    ob = np.argsort(kb, kind='stable'); kbs = kb[ob]
    lo = np.searchsorted(kbs, ka, 'left'); hi = np.searchsorted(kbs, ka, 'right')
    out = []
    for i in np.flatnonzero(hi > lo):
        for j in ob[lo[i]:hi[i]]:
            if np.array_equal(VA[i] + VB[j], target): out.append((int(i), int(j)))
    return out


def group_paf(G, x):
    return np.array([int(np.dot(x, x[G.translate(g)])) for g in range(G.v)], dtype=np.int64)


def verify_pair(G, A, B, T):
    tot = group_paf(G, A) + group_paf(G, B)
    if not (tot[0] == 2 * G.v and np.all(tot[1:] == T)): return False
    if T == 2:
        diff = G.index(G.coords[None, :, :] - G.coords[:, None, :])
        cA, cB = A[diff], B[diff]
        M = np.block([[cA, cB], [-cB.T, cA.T]])
        n = 2 * G.v
        return np.array_equal(M @ M.T, (2 * G.v - 2) * np.eye(n, dtype=M.dtype) + 2 * np.kron(np.eye(2, dtype=M.dtype), np.ones((G.v, G.v), dtype=M.dtype)))
    return True


def freq_table(G, st):
    """S[o, c] = sum over the orbit o of chi_c(x), for a set of nonzero frequencies chi_c (for Z_v one per class of
    <H, -1>, otherwise every nonzero frequency; repeats only cost time in the PSD test)."""
    if len(G.mods) == 1:
        freqs = np.array(st.reps, dtype=np.int64)[:, None]
    else:
        freqs = G.coords[1:]
    ph = (G.coords.astype(np.float64) / G.mods[None, :]) @ freqs.T.astype(np.float64)     # v x F
    E = np.exp(2j * np.pi * ph)
    S = np.zeros((st.O, len(freqs)), dtype=np.complex128)
    np.add.at(S, st.orb, E)
    return S


def psd_survivors(G, st, r, bound, log):
    """Orbit sign vectors with row sum r and |X^(chi)|^2 <= bound for every frequency of freq_table (psd_pool.c)."""
    import subprocess, os, tempfile
    S = freq_table(G, st)
    lines = [f"{st.O} {S.shape[1]} {r} {bound:.6f}", ' '.join(str(int(x)) for x in st.sizes)]
    lines += [' '.join(f"{z.real:.12f} {z.imag:.12f}" for z in S[o]) for o in range(st.O)]
    fd, path = tempfile.mkstemp(suffix='.bin', dir=W.tmpdir()); os.close(fd)
    try:
        out = subprocess.run([W.bin('psd_pool'), path], input='\n'.join(lines) + '\n', capture_output=True, text=True, check=True).stdout
        masks = np.fromfile(path, dtype=np.uint64)
    finally:
        os.remove(path)
    Y = (1 - 2 * ((masks[:, None] >> np.arange(st.O, dtype=np.uint64)) & np.uint64(1))).astype(np.int8)
    return Y


def search(kind, v, name, G, mats, st, maxpool, log, psdmax=0):
    T = {'dopt': 2, 'lp': -2, 'pg': 0}[kind]
    found = []
    for ra, rb in patterns(v, T):
        sizes = [st.pool(ra), st.pool(rb)]
        if min(sizes) == 0 or max(sizes) > max(maxpool, psdmax):
            log(f"  {name} rows ({ra}, {rb}): pools {[int(s) for s in sizes]} skipped"); continue
        if max(sizes) > maxpool and max(sizes) * 0.87 ** st.C > maxpool:   # the PSD test keeps about 0.87 per class
            log(f"  {name} rows ({ra}, {rb}): pools {[int(s) for s in sizes]}, {st.C} classes: too many PSD survivors, skipped"); continue
        t0 = time.time()
        bound = 2 * v - T
        get = lambda r: pool_exact(st, r) if st.pool(r) <= maxpool else psd_survivors(G, st, r, bound, log)
        YA = get(ra); YB = YA if rb == ra else get(rb)
        if max(sizes) > maxpool: log(f"  {name} rows ({ra}, {rb}): PSD test kept {len(YA)} of {int(sizes[0])} and {len(YB)} of {int(sizes[1])}")
        def vals(Y, chunk=1 << 21):
            V = np.empty((len(Y), st.C), dtype=np.int32)
            for c0 in range(0, len(Y), chunk): V[c0:c0 + chunk] = np.rint(st.class_values(Y[c0:c0 + chunk]))
            return V
        VA = vals(YA); VB = VA if rb == ra else vals(YB)
        target = np.full(st.C, T, dtype=np.int32)
        sols = join2(VA, VB, target)
        log(f"  {name} rows ({ra}, {rb}): pools {len(YA)}, {len(YB)}: {len(sols)} solutions [{time.time() - t0:.1f}s]")
        for i, j in sols[:2]:
            A = YA[i].astype(np.int64)[st.orb]; B = YB[j].astype(np.int64)[st.orb]
            ok = verify_pair(G, A, B, T)
            log(f"  hit: {({2: 'D-optimal design', -2: 'Legendre pair', 0: 'periodic Golay pair'})[T]} verified {ok}")
            if ok:
                rec = dict(kind=kind, v=v, group=[int(x) for x in G.mods], structure=name,
                           automorphisms=[np.asarray(M).tolist() for M in mats], rows=[ra, rb],
                           seqs=[''.join('+' if x > 0 else '-' for x in s) for s in (A, B)])
                with open(W.hits('two_circ_hits.jsonl'), 'a') as f: f.write(json.dumps(rec) + '\n')
                found.append(rec)
        if found: break
    return found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('kind', choices=['dopt', 'lp', 'pg']); ap.add_argument('vs')
    ap.add_argument('--maxo', type=int, default=34); ap.add_argument('--maxpool', type=float, default=4e7)
    ap.add_argument('--group', action='store_true', help='non-cyclic groups (split primes) instead of Z_v')
    ap.add_argument('--psdmax', type=float, default=0, help='pools above --maxpool and up to this size go through the PSD test (psd_pool.c)')
    a = ap.parse_args()
    log = lambda s: print(s, flush=True)
    for v in [int(x) for x in a.vs.split(',')]:
        gen = group_structures(v, a.maxo) if a.group else cyclic_structures(v, a.maxo)
        jobs = sorted(gen, key=lambda j: j[3].O)
        log(f"=== {a.kind} v={v}: {len(jobs)} structures with <= {a.maxo} orbits{' (non-cyclic)' if a.group else ''}")
        hit = False
        for name, G, mats, st in jobs:
            try:
                if search(a.kind, v, name, G, mats, st, a.maxpool, log, a.psdmax): hit = True; break
            except MemoryError as e:
                log(f"  {name}: out of memory ({e}); skipped")
        log(f"=== {a.kind} v={v}: {'HIT' if hit else 'none'}")


if __name__ == '__main__':
    main()

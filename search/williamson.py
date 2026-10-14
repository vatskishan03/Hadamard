"""Williamson matrices with multiplier symmetry: four symmetric circulant +-1 matrices A, B, C, D of order n with
A^2 + B^2 + C^2 + D^2 = 4n I (equivalently PAF_A + PAF_B + PAF_C + PAF_D = 4n delta with symmetric sequences).

For every subgroup H of Z_n^* (cyclic subgroups and products of two) the four sequences are taken constant on the orbits
of <H, -1>; for every row pattern r1 <= r2 <= r3 <= r4 (odd, sum of squares 4n) the complete pools are joined exactly
(join4x), up to a cap on the pair work. Hits are rebuilt in the Williamson array and checked (H H^T = 4n I).

  python3 search/williamson.py n1,n2,.. [--maxo 24] [--maxwork 1e11] [--K 1]
Hits go to $HADAMARD_WORK/williamson_hits.jsonl.
"""
import argparse, json, math, os, struct, subprocess, sys, tempfile, time
import numpy as np
from pathlib import Path as _P
sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
import spence_lib as L
import propus as PR
import verify
from skew_ktree import pool_exact


def patterns(n):
    lim = math.isqrt(4 * n); out = []
    for a in range(1, lim + 1, 2):
        for b in range(a, lim + 1, 2):
            for c in range(b, lim + 1, 2):
                d2 = 4 * n - a * a - b * b - c * c
                if d2 < c * c: continue
                d = math.isqrt(d2)
                if d * d == d2 and d % 2: out.append((a, b, c, d))
    return out


def williamson_array(A, B, C, D):
    a, b, c, d = (verify.circulant(x) for x in (A, B, C, D))
    return np.block([[a, b, c, d], [-b, a, -d, c], [-c, d, a, -b], [-d, -c, b, a]])


def search(n, gens, maxwork, K, log):
    G = L.AbGroup([n])
    st = L.Structure(G, [G.perm_from_matrix(np.array([[g % n]])) for g in gens] + [G.neg()])
    pm = n % 4; found = []
    pats = patterns(n)
    log(f"n={n} gens={gens}: orbits {st.O}, classes {st.C}, patterns {len(pats)}")
    cache = {}
    def pool(r):
        if r not in cache:
            Y = pool_exact(st, r)
            V = np.rint((st.class_values(Y) - pm) / 4).astype(np.int16) if len(Y) else np.zeros((0, st.C), np.int16)
            cache[r] = (Y, V)
        return cache[r]
    for q in pats:
        sizes = [st.pool(r) for r in q]
        if min(sizes) == 0: continue
        s = sorted(sizes); work = s[0] * s[3] + s[1] * s[2]
        if work > maxwork:
            log(f"  rows {q}: pools {[int(x) for x in sizes]} work {work:.1e} skipped"); continue
        t0 = time.time()
        lists = [pool(r) for r in q]
        order = sorted(range(4), key=lambda k: sizes[k]); perm = [order[0], order[3], order[1], order[2]]
        fd, path = tempfile.mkstemp(suffix='.bin', dir=W.tmpdir()); os.close(fd)
        try:
            with open(path, 'wb') as f:
                f.write(struct.pack('iiiii', *[len(lists[k][1]) for k in perm], st.C))
                f.write(np.full(st.C, -pm, dtype=np.int16).tobytes())
                for k in perm: f.write(np.ascontiguousarray(lists[k][1]).tobytes())
            out = subprocess.run([W.bin('join4x'), path, str(K)], capture_output=True, text=True, check=True).stdout
        finally:
            os.remove(path)
        total = int(out.strip().split('TOTAL')[-1])
        log(f"  rows {q}: pools {[int(x) for x in sizes]}: {total} solutions [{time.time() - t0:.1f}s]")
        for line in out.splitlines():
            if not line.startswith('SOL'): continue
            ids = list(map(int, line.split()[1:])); idx = [None] * 4
            for slot, k in enumerate(perm): idx[k] = ids[slot]
            seqs = [lists[k][0][idx[k]].astype(np.int64)[st.orb] for k in range(4)]
            tot = verify.paf_sum(seqs)
            ok = tot[0] == 4 * n and not tot[1:].any() and all(np.array_equal(x, x[(-np.arange(n)) % n]) for x in seqs)
            if ok:
                M = williamson_array(*seqs)
                ok = np.array_equal(M @ M.T, 4 * n * np.eye(4 * n, dtype=M.dtype))
            log(f"  hit: Williamson matrices verified {ok}")
            if ok:
                rec = dict(n=n, gens=list(gens), rows=list(q), seqs=[''.join('+' if x > 0 else '-' for x in s_) for s_ in seqs])
                with open(W.hits('williamson_hits.jsonl'), 'a') as fh: fh.write(json.dumps(rec) + '\n')
                found.append(rec)
        if found: break
    return found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('ns'); ap.add_argument('--maxo', type=int, default=24); ap.add_argument('--maxwork', type=float, default=1e11)
    ap.add_argument('--K', type=int, default=1)
    a = ap.parse_args()
    log = lambda s: print(s, flush=True)
    for n in [int(x) for x in a.ns.split(',')]:
        jobs = []
        for S, gs in PR.subgroup_gens(n).items():
            Sm = set(S) | {(-h) % n for h in S}
            O = len({min(t * h % n for h in Sm) for t in range(n)})
            if O <= a.maxo: jobs.append((O, gs))
        jobs.sort()
        log(f"=== n={n}: {len(jobs)} structures with <= {a.maxo} orbits")
        hit = False
        for _, gs in jobs:
            if search(n, gs, a.maxwork, a.K, log): hit = True; break
        log(f"=== n={n}: {'HIT' if hit else 'none'}")


if __name__ == '__main__':
    main()

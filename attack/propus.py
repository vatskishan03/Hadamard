"""Symmetric Hadamard matrices of order 4v from propus-type quadruples (A, B, B, D) on Z_v.

If A and D are symmetric circulant ±1 matrices of order v, B is a circulant ±1 matrix and
A Aᵀ + 2 B Bᵀ + D Dᵀ = 4v I, then with R the back-diagonal identity

        [ A    D    BR   BR ]
    H = [ D   −A    BR  −BR ]
        [ BR   BR  −D   −A  ]
        [ BR  −BR  −A    D  ]

is a symmetric Hadamard matrix of order 4v (every block XR is symmetric for circulant X). In terms of sequences:
PAF_A + 2 PAF_B + PAF_D = 4v δ with A, D symmetric.

The search is exhaustive per multiplier structure: A and D are constant on the orbits of ⟨H, −1⟩ (hence symmetric),
B on the orbits of H (or of ⟨H, −1⟩ with --symB). For every row-sum pattern r_A² + 2 r_B² + r_D² = 4v the complete
pools are listed and joined exactly with join4x on the autocorrelation classes. Hits are rebuilt as matrices and
checked (H = Hᵀ, H Hᵀ = 4v I) before they are recorded in $HADAMARD_WORK/sym_hits.jsonl.

  python3 attack/propus.py v gens [--symB] [--maxpool N] [--K N]      (gens: comma-separated multipliers, or 1)
  python3 attack/propus.py v --sweep [--maxo O]                       (every subgroup structure with <= O orbits)
"""
import argparse, json, math, os, struct, subprocess, sys, tempfile, time
import numpy as np
from pathlib import Path as _P
sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
import spence_lib as L
import count_exact as X
import verify


def pool_vectors(st, r):
    """All sign vectors (int8, one entry per orbit) with row sum r."""
    O = st.O; out = []; step = 1 << 20; sz = st.sizes.astype(np.int64)
    for start in range(0, 1 << O, step):
        idx = np.arange(start, min(start + step, 1 << O), dtype=np.int64)
        Y = (1 - 2 * ((idx[:, None] >> np.arange(O)) & 1)).astype(np.int8)
        out.append(Y[(Y.astype(np.int64) @ sz) == r])
    return np.concatenate(out)


def symmetric_matrix(A, B, D):
    """The propus-type array from first rows A, B, D (A, D symmetric)."""
    cA, cB, cD = (verify.circulant(x) for x in (A, B, D))
    R = verify.back_identity(len(A))
    BR = cB @ R
    return np.block([[cA, cD, BR, BR], [cD, -cA, BR, -BR], [BR, BR, -cD, -cA], [BR, -BR, -cA, cD]])


def symmetric_matrix_v3(A, B, D):
    """Variant needing only A symmetric: [[A, B, B, D], [Bᵀ, DR, −A, −BR], [Bᵀ, −A, −DR, BR], [Dᵀ, −BR, BR, −A]]."""
    cA, cB, cD = (verify.circulant(x) for x in (A, B, D))
    R = verify.back_identity(len(A))
    BR, DR = cB @ R, cD @ R
    return np.block([[cA, cB, cB, cD], [cB.T, DR, -cA, -BR], [cB.T, -cA, -DR, BR], [cD.T, -BR, BR, -cA]])


def patterns(v):
    out = []
    lim = math.isqrt(4 * v)
    for ra in range(1, lim + 1, 2):
        for rb in range(1, lim + 1, 2):
            rest = 4 * v - ra * ra - 2 * rb * rb
            if rest < 1: continue
            rd = math.isqrt(rest)
            if rd * rd == rest and rd % 2: out.append((ra, rb, rd))
    return out


def structures(v, gens, symB):
    G = L.AbGroup([v])
    perm = lambda g: G.perm_from_matrix(np.array([[g % v]]))
    st_ad = L.Structure(G, [perm(g) for g in gens] + [G.neg()])
    st_b = st_ad if symB else L.Structure(G, [perm(g) for g in gens] if gens else [perm(1)])
    assert st_ad.C == st_b.C and np.array_equal(st_ad.cls, st_b.cls), "class partitions differ"
    return st_ad, st_b


def search(v, gens, symB=False, maxpool=3e7, K=1, log=print, freeD=False):
    st_ad, st_b = structures(v, gens, symB)
    st_a, st_d = st_ad, (st_b if freeD else st_ad)
    pm = v % 4
    found = []
    pats = patterns(v)
    log(f"v={v} gens={gens} symB={symB} freeD={freeD}: orbits A {st_a.O}, B {st_b.O}, D {st_d.O}, classes {st_ad.C}, patterns {len(pats)}")
    for ra, rb, rd in pats:
        # signs: negating A, B or D gives an equivalent solution, so positive row sums suffice
        if not freeD and ra > rd: continue                 # A and D interchangeable when both are symmetric
        sizes = [st_a.pool(ra), st_b.pool(rb), st_d.pool(rd)]
        if min(sizes) == 0 or max(sizes) > maxpool:
            log(f"  rows ({ra}, {rb}, {rd}): pools {[int(s) for s in sizes]} skipped"); continue
        t = time.time()
        def vals(st, Y, scale, chunk=1 << 18):
            out = np.empty((len(Y), st.C), dtype=np.int16)
            for c0 in range(0, len(Y), chunk):
                out[c0:c0 + chunk] = scale * np.rint((st.class_values(Y[c0:c0 + chunk].astype(np.float64)) - pm) / 4)
            return out
        Ya = pool_vectors(st_a, ra); Yb = pool_vectors(st_b, rb)
        Yd = Ya if (rd == ra and not freeD) else pool_vectors(st_d, rd)
        Va = vals(st_a, Ya, 1); Vb = vals(st_b, Yb, 2); Vd = Va if (rd == ra and not freeD) else vals(st_d, Yd, 1)
        Z = np.zeros((1, st_ad.C), dtype=np.int16)
        lists = [Va, Vd, Vb, Z]                     # join4x stores the smaller pair product: A x D against B
        fd, path = tempfile.mkstemp(suffix='.bin', dir=W.tmpdir()); os.close(fd)
        try:
            with open(path, 'wb') as f:
                f.write(struct.pack('iiiii', *[len(x) for x in lists], st_ad.C))
                f.write(np.full(st_ad.C, -pm, dtype=np.int16).tobytes())
                for x in lists: f.write(np.ascontiguousarray(x).tobytes())
            out = subprocess.run([W.bin('join4x'), path, str(K)], capture_output=True, text=True, check=True).stdout
        finally:
            os.remove(path)
        total = int(out.strip().split('TOTAL')[-1])
        log(f"  rows ({ra}, {rb}, {rd}): pools {[int(s) for s in sizes]}: {total} solutions [{time.time() - t:.1f}s]")
        for line in out.splitlines():
            if not line.startswith('SOL'): continue
            ia, idd, ib, _ = map(int, line.split()[1:])
            A = Ya[ia].astype(np.int64)[st_a.orb]; B = Yb[ib].astype(np.int64)[st_b.orb]; D = Yd[idd].astype(np.int64)[st_d.orb]
            tot = verify.paf_sum([A, B, B, D])
            if not (tot[0] == 4 * v and not tot[1:].any()): log('  WRONG autocorrelations'); continue
            Hm = symmetric_matrix_v3(A, B, D) if freeD else symmetric_matrix(A, B, D)
            ok = np.array_equal(Hm, Hm.T) and np.array_equal(Hm @ Hm.T, 4 * v * np.eye(4 * v, dtype=Hm.dtype))
            log(f"  hit: symmetric Hadamard {ok}")
            if ok:
                rec = dict(n=v, gens=list(gens), symB=symB, freeD=freeD, rows=[int(A.sum()), int(B.sum()), int(D.sum())],
                           seqs=[''.join('+' if x > 0 else '-' for x in s) for s in (A, B, D)])
                with open(W.hits('sym_hits.jsonl'), 'a') as fh: fh.write(json.dumps(rec) + '\n')
                found.append(rec)
        if found: break
    return found


def subgroup_gens(v):
    """One generating set per subgroup of Z_v^* (cyclic subgroups and products of two)."""
    units = [u for u in range(1, v) if math.gcd(u, v) == 1]
    def gen(gs):
        Hs = {1}
        while True:
            H2 = Hs | {h * g % v for h in Hs for g in gs}
            if H2 == Hs: return frozenset(Hs)
            Hs = H2
    subs = {}
    cyc = {}
    for u in units:
        S = gen([u])
        if S not in cyc: cyc[S] = [u]
    for S, g in cyc.items(): subs.setdefault(S, g)
    keys = list(cyc.items())
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            S = gen(keys[i][1] + keys[j][1])
            if S not in subs: subs[S] = keys[i][1] + keys[j][1]
    return subs


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('v', type=int); ap.add_argument('gens', nargs='?', default='1')
    ap.add_argument('--symB', action='store_true'); ap.add_argument('--maxpool', type=float, default=3e7)
    ap.add_argument('--K', type=int, default=1); ap.add_argument('--sweep', action='store_true')
    ap.add_argument('--maxo', type=int, default=30); ap.add_argument('--mino', type=int, default=1)
    ap.add_argument('--freeD', action='store_true', help='only A symmetric (array of symmetric_matrix_v3)')
    a = ap.parse_args()
    log = lambda s: print(s, flush=True)
    if a.sweep:
        v = a.v
        jobs = []
        for S, gs in subgroup_gens(v).items():
            Hm = set(S) | {(-h) % v for h in S}
            Ob = len({min(t * h % v for h in S) for t in range(v)})
            Oad = len({min(t * h % v for h in Hm) for t in range(v)})
            if Oad < a.mino or max(Oad, Ob) > a.maxo: continue
            jobs.append((Ob + Oad, gs))
        jobs.sort()
        log(f"v={v}: {len(jobs)} structures")
        for _, gs in jobs:
            for symB, freeD in ((True, False), (False, False), (False, True)):
                try:
                    if search(v, gs, symB, a.maxpool, a.K, log, freeD): log(f"HIT v={v} gens={gs} symB={symB} freeD={freeD}"); sys.exit(0)
                except AssertionError as e:
                    log(f"  skipped: {e}")
    else:
        gens = [int(x) for x in a.gens.split(',')]
        search(a.v, gens, a.symB, a.maxpool, a.K, log, a.freeD)

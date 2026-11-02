"""Cross-check of p3.c against an independent reference solver on many small structures.

Reference (pure numpy/dict, no C, no hashing tricks): all A, B, D sign vectors with the row sums, class values from the
S-tensors (themselves verified against direct PAFs in p3struct.selfcheck), exact dictionary join on the whole class vector.
For every structure and row pattern the C engine must return
   * unreduced mode (no canonical filters; PSD filters on or off): exactly the reference solution set;
   * reduced mode: exactly the reference solutions whose B is canonical under Mult and D under reversal of D;
   * for every ns in {1,2,3}, threads in {1,3}, and shards {1, 3}: the same set (union over shards).
Every reference solution is also re-checked by direct periodic autocorrelations of the sequences on Z_v.

  python3.11 p3test.py [--maxpool 40000] [--vmax 100] [--quick]
"""
import argparse, itertools, math, sys, tempfile, time
from math import gcd
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from p3struct import Structure, subgroup, paf
import p3run as P


def all_vectors(blk, r, cap):
    """All sign vectors (n, O) with row sum r (enumerated by bit masks; O <= 24 here)."""
    O = blk.O
    n = 1 << O
    if n > cap:
        return None
    m = np.arange(n, dtype=np.int64)
    bits = ((m[:, None] >> np.arange(O)[None, :]) & 1).astype(np.int64)
    Y = 2 * bits - 1
    s = Y @ blk.sizes
    sel = np.flatnonzero(s == r)
    return m[sel], Y[sel].astype(np.int8)


def reference_solutions(st, rows, tgt=None):
    ra, rb, rd = rows
    tgt = np.zeros(st.C, dtype=np.int64) if tgt is None else np.asarray(tgt)
    A = all_vectors(st.A, ra, 1 << 22)
    Bv = all_vectors(st.B, rb, 1 << 22)
    D = all_vectors(st.B, rd, 1 << 22)
    if A is None or Bv is None or D is None or len(A[0]) == 0 or len(Bv[0]) == 0 or len(D[0]) == 0:
        return None
    FA, FB, FD = st.A.values(A[1]), st.B.values(Bv[1]), st.B.values(D[1])
    w = np.random.default_rng(7).integers(1, 1 << 40, size=st.C, dtype=np.int64)
    Dh = FD @ w
    order = np.argsort(Dh, kind='stable')
    Dhs = Dh[order]
    Bh = 2 * (FB @ w)
    sols = set()
    for i in range(len(A[0])):
        nh = (tgt - FA[i]) @ w - Bh
        pos = np.searchsorted(Dhs, nh)
        pos = np.minimum(pos, len(Dhs) - 1)
        cand = np.flatnonzero(Dhs[pos] == nh)
        for j in cand:
            p = pos[j]
            while p < len(Dhs) and Dhs[p] == nh[j]:
                d = order[p]
                if np.array_equal(FA[i] + 2 * FB[j] + FD[d], tgt):
                    sols.add((int(A[0][i]), int(Bv[0][j]), int(D[0][d])))
                p += 1
    return sols


def structures(vmax, maxpool):
    out = []
    for v in range(15, vmax + 1, 2):
        units = [t for t in range(2, v) if gcd(t, v) == 1]
        seen = set()
        for u in units:
            H = tuple(subgroup(v, [u]))
            if H in seen or len(H) == 1 or (v - 1) in H:
                continue
            seen.add(H)
            try:
                st = Structure(v, [u])
            except AssertionError:
                continue
            if st.A.O > 22 or st.B.O > 22:
                continue
            for p in st.feasible_patterns():
                if max(p['pools'][1:]) <= maxpool and p['pools'][0] * p['pools'][1] <= 4e6:
                    out.append((v, u, p['rows']))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--maxpool', type=int, default=40000)
    ap.add_argument('--vmax', type=int, default=100)
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--limit', type=int, default=10 ** 9)
    a = ap.parse_args()
    binary = P.build_binary()
    cases = structures(a.vmax, a.maxpool)[:a.limit]
    print(f'{len(cases)} (v, H, rows) cases', flush=True)
    tmp = Path(tempfile.mkdtemp(prefix='p3test_'))
    rng = np.random.default_rng(5)
    nsol_total = 0
    nfail = 0
    cache = {}
    t0 = time.process_time()
    for ci, (v, u, rows) in enumerate(cases):
        if (v, u) not in cache:
            cache[(v, u)] = Structure(v, [u])
        st = cache[(v, u)]
        ref = reference_solutions(st, rows)
        if ref is None:
            continue
        # direct re-check of every reference solution
        for (ma, mb, md) in list(ref)[:50]:
            assert P.exact_check(st, *P.masks_to_seqs(st, ma, mb, md)), 'reference solution fails the direct PAF check'
        refc = {s for s in ref if P.is_canonical(st, *s)}
        # losslessness of the reduction, independent of the C engine: the Gamma-orbits of the canonical solutions cover all
        expanded = set()
        for s0 in refc:
            expanded |= P.orbit_group_images(st, *s0)
        if expanded != ref:
            nfail += 1
            print(f'FAIL v={v} H={len(st.H)} rows={rows}: Gamma-expansion of canonical solutions ({len(expanded)}) != all ({len(ref)})', flush=True)
        nsol_total += len(ref)
        variants = [(False, (True, True, True)), (False, (False, False, False)), (True, (True, True, True))]
        if a.quick:
            variants = variants[:1] + variants[2:]
        for reduce, psd in variants:
            pf = tmp / 'prob.txt'
            P.write_problem(st, rows, pf, canonB=reduce, canonD=reduce, psd=psd, seed=int(rng.integers(1, 1 << 30)))
            want = refc if reduce else ref
            for ns, thr, shards in ((1, 1, 1), (2, 3, 1), (3, 2, 3)):
                got = set()
                for sh in range(shards):
                    res = P.run_p3(pf, binary, threads=thr, ns=ns, shard=(sh, shards) if shards > 1 else None)
                    got |= set(res['hits'])
                    assert res['done'] is not None and int(res['done']['badmatches']) == 0, 'bad hash matches (should be impossible)'
                if got != want:
                    nfail += 1
                    print(f'FAIL v={v} H={len(st.H)} rows={rows} reduce={reduce} psd={psd} ns={ns} thr={thr} shards={shards}: '
                          f'C {len(got)} vs ref {len(want)}', flush=True)
        if len(ref) and ci % 5 == 0:
            print(f'  v={v:3d} |H|={len(st.H):2d} OA={st.A.O:2d} OB={st.B.O:2d} rows={rows}: ref {len(ref)} solutions, canonical {len(refc)}', flush=True)
    print(f'total reference solutions over all cases: {nsol_total}; failures: {nfail}; python cpu {time.process_time() - t0:.1f}s')
    sys.exit(1 if nfail else 0)


if __name__ == '__main__':
    main()

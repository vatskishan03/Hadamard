"""Validation of the three-list engine on the repository's two known propus matrices (orders 532 and 876).

For each of symmetric/532 (v = 133, H = <26>, rows (11, 5, 19); 72 solutions = 2 orbits of Gamma, NOTES.md section 30) and
symmetric/876 (v = 219, H = <4>, rows (3, 17, 17); 288 solutions = 1 orbit):
  1. unreduced run (no canonical filters, PSD filters on): must return all solutions, including the repository's;
  2. every hit is re-checked by direct periodic autocorrelations of the sequences on Z_v;
  3. reduced run (canonical B and D): must return exactly the canonical members, and their Gamma-orbits must cover the
     unreduced set;
  4. unreduced run with PSD filters off: same set (the filters lose nothing).
The repository is only read (certificate.json).

  python3.11 p3validate.py [--threads 4] [--only 532|876]
"""
import argparse, json, sys, tempfile, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from p3struct import Structure
import p3run as P

REPO = HERE.parent.parent          # the repository root when this directory is search/propus3


def repo_masks(st, cert):
    v = st.v
    seqs = []
    for reps in cert['blocks']:
        x = np.ones(v, dtype=np.int64)
        for r in reps:
            for h in st.H:
                x[r * h % v] = -1
        seqs.append(x)
    assert P.exact_check(st, *seqs), 'repository certificate fails the direct check?!'
    return P.seqs_to_masks(st, *seqs), seqs


def one(order, v, gens, rows, expect_all, expect_reduced, threads, tmp, psd_off_check=True):
    print(f'=== order {order} (v={v}, H={gens}, rows={rows})', flush=True)
    st = Structure(v, gens)
    print(f'  |H|={len(st.H)} OA={st.A.O} OB={st.B.O} C={st.C} |Mult|={len(st.mult_reps)} pools={[st.A.pool(rows[0]), st.B.pool(rows[1]), st.B.pool(rows[2])]}', flush=True)
    cert = json.loads((REPO / 'symmetric' / str(order) / 'certificate.json').read_text())
    assert cert['v'] == v and cert['rowsums'] == list(rows)
    rm, rseqs = repo_masks(st, cert)
    print(f'  repository solution as orbit masks: {tuple(hex(x) for x in rm)}', flush=True)
    binary = P.build_binary()
    pf = tmp / f'p{order}.txt'
    out = {}
    variants = [('unreduced+psd', False, (True, True, True)), ('reduced+psd', True, (True, True, True))]
    if psd_off_check:
        variants.append(('unreduced, psd off', False, (False, False, False)))
    for label, reduce, psd in variants:
        P.write_problem(st, rows, pf, canonB=reduce, canonD=reduce, psd=psd, seed=11)
        t = time.time()
        res = P.run_p3(pf, binary, threads=threads, ns=2)
        hits = set(res['hits'])
        out[label] = hits
        print(f'  [{label}] lists A,B,D = {res["lists"]}, hits {len(hits)}, pairs {res["done"]["pairs"]}, matches {res["done"]["matches"]}, '
              f'badmatches {res["done"]["badmatches"]}, wall {time.time() - t:.1f}s', flush=True)
        assert int(res['done']['badmatches']) == 0
        for h in hits:
            assert P.exact_check(st, *P.masks_to_seqs(st, *h)), 'engine hit fails the direct check'
    full = out['unreduced+psd']
    ok = True
    ok &= len(full) == expect_all
    print(f'  unreduced count {len(full)} (expected {expect_all}): {"ok" if len(full) == expect_all else "MISMATCH"}')
    if psd_off_check:
        ok &= out['unreduced, psd off'] == full
        print(f'  PSD filters lossless (same set with filters off): {out["unreduced, psd off"] == full}')
    ok &= rm in full
    print(f'  repository solution found: {rm in full}')
    canon = {h for h in full if P.is_canonical(st, *h)}
    red = out['reduced+psd']
    ok &= red == canon
    print(f'  reduced run == canonical members of the full set: {red == canon}  (count {len(red)}, expected {expect_reduced})')
    ok &= len(red) == expect_reduced
    exp = set()
    for h in red:
        exp |= P.orbit_group_images(st, *h)
    ok &= exp == full
    print(f'  Gamma-expansion of the reduced hits == full set: {exp == full}; |Gamma| generic = {P.gamma_order(st)}')
    # the repository solution's own orbit
    img = P.orbit_group_images(st, *rm)
    print(f'  Gamma-orbit of the repository solution has {len(img)} members; subset of hits: {img <= full}')
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--only', type=int, default=0)
    a = ap.parse_args()
    tmp = Path(tempfile.mkdtemp(prefix='p3val_'))
    ok = True
    if a.only in (0, 532):
        ok &= one(532, 133, [26], (11, 5, 19), 72, 2, a.threads, tmp)
    if a.only in (0, 876):
        ok &= one(876, 219, [4], (3, 17, 17), 288, 9, a.threads, tmp, psd_off_check=False)
    print('ALL OK' if ok else 'FAILED')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()

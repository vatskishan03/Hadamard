"""test_plain.py -- validation of the corrected solver plain_spence3.solve (search/NOTES.md section 92) on small odd v:
 (a) the solver as shipped, (b) the solver with forced bucket splitting (maxpairs = 77: nparts > 1), (c) the same with forced chunked pair generation (CH = 40 in place of 1,500,000: the source text is edited at run time,
 so the shipped file is not changed), and (d) an independent count: PAF vectors of the +1-row-sum sequences by direct rolling on Z_v, ordered pairs in a dictionary, no hashing, no buckets (only for lists L <= 400).
All four must give the same number of ordered quadruples (A, B, C, D) of H-invariant +-1 sequences with all row sums +1 and PAF = -4 off 0 (the printed 'count_rho_plus' of plain_spence3.py).
usage: test_plain.py [vmax [Omax]]    (default 75, 11)    prints one line per shape with a nonzero count and a summary"""
import sys, os, time, collections, types
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import plain_spence3 as ps
src = open(os.path.join(HERE, 'plain_spence3.py')).read()
assert src.count('CH = 1500000') == 1
small = types.ModuleType('plain_spence3_ch40'); exec(compile(src.replace('CH = 1500000', 'CH = 40'), 'plain_spence3_ch40', 'exec'), small.__dict__)

def solver_count(mod, v, H, maxpairs=None):
    lab, reps, O, size, cls = mod.structure(v, H)
    Ms = mod.forms(v, lab, reps, O, cls); X = mod.enum_list(O, size)
    if len(X) == 0: return 0, 0
    Xi = X.astype(np.int64)
    V = np.stack([np.einsum('ni,ij,nj->n', Xi, M, Xi) for M in Ms], 1)
    Vu, inv = np.unique(V, axis=0, return_inverse=True); cnt = np.bincount(inv.ravel())
    tgt = np.full(len(cls), -4, dtype=np.int64)
    tot, found = mod.solve(Vu, cnt, tgt) if maxpairs is None else mod.solve(Vu, cnt, tgt, maxpairs=maxpairs)
    return int(tot), len(X)

def indep_count(v, H):
    lab, reps, O, size, cls = ps.structure(v, H)
    X = ps.enum_list(O, size); S = X[:, lab].astype(np.int64)
    assert (S.sum(1) == 1).all()
    V = np.stack([(S * np.roll(S, -g, axis=1)).sum(1) for g in range(1, v)], 1)
    P = collections.Counter(map(bytes, (V[:, None, :] + V[None, :, :]).reshape(-1, v - 1).astype(np.int16)))
    T = np.full(v - 1, -4, dtype=np.int16)
    tot = 0
    for k, c in P.items():
        a = np.frombuffer(k, dtype=np.int16)
        tot += c * P.get(bytes(T - a), 0)
    return tot

vmax = int(sys.argv[1]) if len(sys.argv) > 1 else 75; Omax = int(sys.argv[2]) if len(sys.argv) > 2 else 11
nshape = nzero = nind = bad = nnz = 0; t0 = time.time()
for v in range(9, vmax + 1, 2):
    for H in ps.all_subgroups(v):
        lab, reps, O, size, cls = ps.structure(v, H)
        if not (2 <= O <= Omax): continue
        a, L = solver_count(ps, v, H)
        if L == 0: continue
        b, _ = solver_count(ps, v, H, maxpairs=77); c, _ = solver_count(small, v, H, maxpairs=77)
        res = [a, b, c]
        if L <= 400: res.append(indep_count(v, H)); nind += 1
        nshape += 1; nnz += (a > 0)
        if len(set(res)) != 1: bad += 1; print('MISMATCH v', v, 'h', len(H), 'O', O, 'L', L, res, flush=True)
        elif a: print('v %d |H| %d O %d L %d: count %d (all %d methods agree)' % (v, len(H), O, L, a, len(res)), flush=True)
print('test_plain: %d shapes with L > 0 (v = 9..%d, 2 <= O <= %d), %d with a nonzero count, %d with the independent count, %d mismatches  [%.0f s]' % (nshape, vmax, Omax, nnz, nind, bad, time.time() - t0))

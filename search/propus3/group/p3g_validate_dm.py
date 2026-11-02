"""Validation of the dmult reduction (D minimal under Mult x reversal, B <= reversal) against brute force."""
import sys, os, itertools, tempfile
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p3g
from p3g import GStruct
import importlib.util
spec = importlib.util.spec_from_file_location('pv', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'p3g_validate.py'))
# reuse brute() from the validator without running it: copy minimal brute force here
import propus_group as PG
work = tempfile.mkdtemp(prefix='p3gdm_')

def all_vectors(blk, r):
    O = blk.O; sizes = [int(x) for x in blk.sizes]; outs = []
    def rec(j, cur, s):
        if j == O:
            if s == r: outs.append(list(cur))
            return
        if abs(r - s) > sum(sizes[j:]): return
        for y in (1, -1):
            cur.append(y); rec(j + 1, cur, s + y * sizes[j]); cur.pop()
    rec(0, [], 0)
    return np.array(outs, dtype=np.int8)

def mask(y): return sum(1 << i for i in range(len(y)) if y[i] > 0)

def brute(st, rows):
    G = st.G; ra, rb, rd = rows
    YA = all_vectors(st.A, ra); YB = all_vectors(st.B, rb); YD = all_vectors(st.B, rd)
    pafs = lambda blk, Y: np.array([PG.group_paf(G, Y[i][blk.orb].astype(np.int64)) for i in range(len(Y))], dtype=np.int64)
    PA, PB, PD = pafs(st.A, YA), pafs(st.B, YB), pafs(st.B, YD)
    dd = {}
    for j in range(len(YD)): dd.setdefault(PD[j, 1:].tobytes(), []).append(j)
    sols = set()
    for i in range(len(YA)):
        for k in range(len(YB)):
            js = dd.get((-PA[i, 1:] - 2 * PB[k, 1:]).astype(np.int64).tobytes())
            if js:
                for j in js: sols.add((mask(YA[i]), mask(YB[k]), mask(YD[j])))
    return sols

tests = [(25, [(0, 1)]), (25, [(8, 0)]), (27, [(0, 1)]), (49, [(2, 1)]), (49, [(16, 0)]), (63, [(0, 2, 1, 0)]), (45, [(0, 1)]), (75, [(2, 2, 0, 0)]), (25, [(0, 1), (8, 0)])]
ok_all = True; ncase = 0; nsol = 0
for v, gens in tests:
    try:
        st = GStruct(v, hgens=gens)
    except Exception as ex:
        print('skip', v, gens, ex); continue
    for pat in st.feasible_patterns():
        rows = pat['rows']
        if pat['pools'][0] * pat['pools'][1] > 6e5 or pat['pools'][2] > 3e5: continue
        truth = brute(st, rows)
        canon = {s for s in truth if p3g.is_canonical_dmult(st, *s)}
        r1 = set(p3g.run_dmult(st, rows, work, threads=1, ns=2, tag='a')['hits'])
        r2 = set(p3g.run_dmult(st, rows, work, threads=3, ns=3, tag='b')['hits'])
        r3 = set(p3g.run_dmult(st, rows, work, threads=1, ns=2, nopsd=True, tag='c')['hits'])
        sh = set()
        for i in range(3): sh |= set(p3g.run_dmult(st, rows, work, threads=2, ns=2, shard=(i, 3), tag='d%d' % i)['hits'])
        exp = set()
        for s in r1: exp |= st.gamma_images(*s)
        ok = (r1 == canon and r2 == canon and r3 == canon and sh == canon and exp == truth)
        for h in r1:
            A, B, D = st.masks_to_seqs(*h)
            ok = ok and st.exact_check(A, B, D) and st.matrix_check(A, B, D)
        ncase += 1; nsol += len(truth); ok_all &= ok
        # compare with the default reduction sizes
        n_def = len({s for s in truth if st.is_canonical(*s)})
        print(f'v={v} gens={gens} |H|={st.Hord} |mult|={len(st.mult_reps)} fix={len(st.fix)} rows={rows}: brute {len(truth)}, default-canonical {n_def}, dmult-canonical {len(canon)} | p3 dmult {len(r1)} 3thr {len(r2)} nopsd {len(r3)} shards {len(sh)} gamma-exp {len(exp)} -> {"OK" if ok else "FAIL"}', flush=True)
print(f'cases {ncase}, solutions {nsol}, ALL OK' if ok_all else 'FAILED'); sys.exit(0 if ok_all else 1)

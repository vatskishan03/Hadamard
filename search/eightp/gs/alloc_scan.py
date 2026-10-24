"""List the defect allocations (n0<=n1<=n2<=n3, optionally m-splits) of a shape that pass the exact identity AND the exact row-sum pattern filter (gsk.row_patterns).
usage: alloc_scan.py MODE Kmods r hg T Emax [mfix]   ->  lines  'n0:n1:n2:n3  E  npatterns  log10(prod of block spaces)'
(with mfix = 1 the m-splits (m0..m3) with sum = e_m are enumerated per allocation as well: 'n.m' tokens, usable directly as gs_samp.py allocs)"""
import sys, os, json, itertools, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gs_samp import *
mode = sys.argv[1]; Kmods = json.loads(sys.argv[2]); r = int(sys.argv[3]); hg = json.loads(sys.argv[4]); T = json.loads(sys.argv[5]); Emax = int(sys.argv[6])
mfix = len(sys.argv) > 7 and sys.argv[7] == '1'
sh = make_shape(mode, Kmods, r, hg, T)
W = [n for n in sh.kCD.weights() if n <= Emax]
ds = [d for d in sh.feasible_totals() if d <= Emax]
print('# shape', mode, Kmods, r, hg, T, 'N', sh.N, 'nu', sh.nu, 'weights', sh.kCD.weights()[:12], 'feasible E', ds, 'T2', sh.T2, 'Rsq', sh.Rsq, flush=True)
for al in itertools.combinations_with_replacement(W, 4):
    E = sum(al)
    if E not in ds: continue
    em = ((sh.r + 1) * E - 4 * sh.N - int(sh.phi[0])) // 4
    t0 = time.time()
    if not mfix:
        P = row_patterns(sh, al, (-1,) * 4)
        sp = sum(math.log10(max(space_size(sh, n), 1)) for n in al)
        print(':'.join(map(str, al)), E, 'e_m', em, 'patterns', None if P is None else len(P), f'log10space {sp:.1f}', f'{time.time()-t0:.0f}s', flush=True)
    else:
        seen = set()
        for ms in itertools.product(*[range(0, n + 1) for n in al]):
            if sum(ms) != em: continue
            key = tuple(sorted(zip(al, ms)))
            if key in seen: continue
            seen.add(key)
            P = row_patterns(sh, al, ms)
            if P is None or len(P):
                print(':'.join(f'{n}.{m}' for n, m in zip(al, ms)), 'E', E, 'patterns', None if P is None else len(P), flush=True)

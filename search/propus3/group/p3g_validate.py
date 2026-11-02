"""Validation of p3g against an independent brute force (direct group PAF + dictionary join) on small structures.
For every (v, e) in the scan with feasible patterns and small pools, and every pattern:
  1. class tensors / theta (PSD) vs direct group PAF and direct DFT on random sign vectors  (selfcheck)
  2. all solutions by brute force = dictionary join on full group PAF vectors of the explicit pools (no p3, no repo join)
  3. p3 --noreduce --nopsd == brute force (as sets of orbit-mask triples); p3 --noreduce (PSD on) == same;
     p3 reduced (canonical B under Aut'/H, canonical D under reversal, PSD on) == canonical members of the brute-force set;
     Gamma-expansion of the reduced hits == brute-force set
  4. threads/ns/shards: 3 threads ns=2 and 3 shards give the same hit set
usage: python3 p3g_validate.py [--vs 25,27,...] [--maxpairs 3e6]"""
import sys, os, json, itertools, time, argparse, tempfile
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p3g
from p3g import GStruct
import propus_group as PG
import p3run as P

ap = argparse.ArgumentParser()
ap.add_argument('--vs', default='25,27,45,49,63,75,81,99')
ap.add_argument('--maxpairs', type=float, default=3e6); ap.add_argument('--cmax', type=int, default=14)
ap.add_argument('--maxsol', type=int, default=100000); ap.add_argument('--mode', default='both')
a = ap.parse_args()
rng = np.random.default_rng(1)
work = tempfile.mkdtemp(prefix='p3gval_')


def all_vectors(blk, r, cap):
    """all sign vectors (n, O) with sum_i size_i y_i = r (meet in the middle over two halves is overkill: DFS bit loop)."""
    O = blk.O; sizes = [int(x) for x in blk.sizes]
    outs = []
    def rec(j, cur, s):
        if j == O:
            if s == r: outs.append(list(cur))
            return
        rest = sum(sizes[j:])
        if abs(r - s) > rest: return
        for y in (1, -1):
            cur.append(y); rec(j + 1, cur, s + y * sizes[j]); cur.pop()
    rec(0, [], 0)
    assert len(outs) <= cap
    return np.array(outs, dtype=np.int8)


def mask(y):
    return sum(1 << i for i in range(len(y)) if y[i] > 0)


def selfcheck(st):
    G, v = st.G, st.v
    for blk, name in ((st.A, 'A'), (st.B, 'B')):
        Y = (rng.integers(0, 2, size=(20, blk.O)) * 2 - 1).astype(np.int8)
        for y in Y:
            x = blk.seq(y).astype(np.int64)
            paf = PG.group_paf(G, x)
            # class values from tensors: F_c = const0 + sum_{i<j} T[c][i][j] y_i y_j
            F = blk.const0 + 0.5 * np.einsum('cij,i,j->c', blk.T, y.astype(np.float64), y.astype(np.float64))
            assert np.array_equal(np.rint(F).astype(np.int64), paf[st.creps]), f'{name}: class values differ from direct PAF'
            for c in range(st.C):
                cols = np.flatnonzero(st.stad.cls == c)
                assert np.all(paf[cols] == paf[cols[0]]), f'{name}: PAF not class-constant'
            # PSD: direct DFT at the character reps vs theta
            coords = G.coords; mods = G.mods.astype(np.float64)
            dft = np.zeros(v, dtype=np.complex128)
            Z = np.array([np.exp(2j * np.pi * ((coords * coords[g]) / mods).sum(1)) for g in range(v)])
            ps = np.abs(Z @ x) ** 2
            # check theta against the direct values at ALL characters: each character value must equal v + theta[c'] F for its class
            Pm = v + st.theta @ paf[st.creps]
            vals = sorted(np.round(Pm, 6).tolist()); allv = ps
            # every theta-predicted value must occur among direct values with right multiplicity classes
            for pv in Pm:
                assert np.min(np.abs(ps - pv)) < 1e-6, f'{name}: theta value {pv} not a direct PSD value'
    return True


def brute(st, rows):
    G, v = st.G, st.v
    ra, rb, rd = rows
    YA = all_vectors(st.A, ra, 1 << 22); YB = all_vectors(st.B, rb, 1 << 22); YD = all_vectors(st.B, rd, 1 << 22)
    def pafs(blk, Y):
        X = Y[:, blk.orb].astype(np.int64)
        return np.array([PG.group_paf(G, x) for x in X], dtype=np.int64)
    PA, PB, PD = pafs(st.A, YA), pafs(st.B, YB), pafs(st.B, YD)
    dd = {}
    for j in range(len(YD)):
        dd.setdefault(PD[j, 1:].tobytes(), []).append(j)
    sols = set()
    for i in range(len(YA)):
        for k in range(len(YB)):
            tgt = (-PA[i, 1:] - 2 * PB[k, 1:])
            js = dd.get(tgt.astype(np.int64).tobytes())
            if js:
                for j in js:
                    sols.add((mask(YA[i]), mask(YB[k]), mask(YD[j])))
    return sols, (len(YA), len(YB), len(YD))


nok = 0; nfail = 0; tot_cases = 0; tot_sols = 0


def test_struct(st, tag):
    global nok, nfail, tot_cases, tot_sols
    v = st.v
    feas = st.feasible_patterns()
    feas = [f for f in feas if f['pools'][0] * f['pools'][1] <= a.maxpairs and f['pools'][2] <= 3e5]
    if not feas: return
    selfcheck(st)
    for f in feas:
        rows = f['rows']
        truth, sizes = brute(st, rows)
        for s_ in truth:
            A, B, D = st.masks_to_seqs(*s_)
            assert st.exact_check(A, B, D), 'brute-force solution fails exact check'
        tot_cases += 1; tot_sols += len(truth)
        r0 = p3g.run(st, rows, work, threads=1, ns=2, noreduce=True, nopsd=True)
        r1 = p3g.run(st, rows, work, threads=1, ns=2, noreduce=True, nopsd=False)
        r2 = p3g.run(st, rows, work, threads=1, ns=2, noreduce=False, nopsd=False)
        r3 = p3g.run(st, rows, work, threads=3, ns=3, noreduce=False, nopsd=False)
        s0, s1, s2, s3 = (set(r['hits']) for r in (r0, r1, r2, r3))
        canon = {s_ for s_ in truth if st.is_canonical(*s_)}
        exp = set()
        for s_ in s2: exp |= st.gamma_images(*s_)
        ok = (s0 == truth and s1 == truth and s2 == canon and s3 == canon and exp == truth)
        sh = set()
        for i in range(3):
            sh |= set(p3g.run(st, rows, work, threads=2, ns=2, noreduce=False, nopsd=False, shard=(i, 3))['hits'])
        ok = ok and sh == canon
        for h in s2:
            A, B, D = st.masks_to_seqs(*h)
            ok = ok and st.exact_check(A, B, D) and st.matrix_check(A, B, D)
        nok += ok; nfail += (not ok)
        print(f"v={v} {tag} |H|={st.Hord} OA={st.A.O} OB={st.B.O} C={st.C} |mult|={len(st.mult_reps)} fix={len(st.fix)} rows={rows} pools={sizes}: "
              f"brute {len(truth)} canon {len(canon)} | p3 nr/np {len(s0)} nr {len(s1)} red {len(s2)} red3t {len(s3)} shards {len(sh)} gamma-exp {len(exp)} -> {'OK' if ok else 'FAIL'}", flush=True)


for v in [int(x) for x in a.vs.split(',')]:
    G, parts = p3g.build_group(v)
    ns = [p ** aa - 1 for p, aa, M in parts]
    seen = set()
    if a.mode in ('cyclic', 'both'):
        for e in itertools.product(*[range(n) for n in ns]):
            if not any(e): continue
            try:
                st = GStruct(v, e)
            except AssertionError:
                continue
            key = (tuple(sorted(np.bincount(st.stb.orb).tolist())), tuple(sorted(np.bincount(st.stad.orb).tolist())), len(st.mult_reps))
            if key in seen or st.C > a.cmax: continue
            seen.add(key)
            test_struct(st, f'cyc e={list(e)}')
    if a.mode in ('semi', 'both'):
        amb = p3g.Ambient(G, parts); els = amb.elements(); neg = amb.neg_elem()
        cyc = {}
        for g in els:
            H = frozenset(amb.closure([g], limit=1 << 14))
            if neg in H: continue
            cyc.setdefault(H, g)
        subs = dict((H, [g]) for H, g in cyc.items())
        cl = list(cyc.items())
        for i, (H1, g1) in enumerate(cl):
            for H2, g2 in cl[i + 1:]:
                H = amb.closure([g1, g2], limit=1 << 12)
                if H is None or neg in H: continue
                H = frozenset(H)
                if H not in subs: subs[H] = [g1, g2]
        for H, gens in subs.items():
            if len(H) == 1: continue
            try:
                st = GStruct(v, hgens=gens)
            except AssertionError:
                continue
            key = (len(H), tuple(sorted(np.bincount(st.stb.orb).tolist())), tuple(sorted(np.bincount(st.stad.orb).tolist())), len(st.mult_reps), len(st.fix))
            if key in seen or st.C > a.cmax: continue
            seen.add(key)
            test_struct(st, f'semi gens={gens}')
print(f'cases {tot_cases}, solutions in total {tot_sols}, passed {nok}, failed {nfail}')
sys.exit(0 if nfail == 0 else 1)

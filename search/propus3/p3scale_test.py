"""Production-configuration checks at the dimensions of v = 239 (O = 35 orbits, C = 17 classes, CP = 32, 5-byte masks),
which the small-structure suite (C <= 16 mostly, O <= 27) does not exercise:
  1. list filters: random B and D masks from the true pools are run through the engine's filters (override lists, pass-1
     counts) and compared with numpy: PSD bound, canonical form under Mult (B) / reversal (D), both together;
  2. DFS enumeration: the engine's leaf count over a slice of the DFS tasks equals the number of completions of the same
     prefixes computed from the row-sum dynamic program.
   python3.11 p3scale_test.py"""
import subprocess, sys, tempfile, time
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import p3run as P
from p3struct import Structure

st = Structure(239, [10]); rows = (15, 1, 27); v = 239
rng = np.random.default_rng(12)
tmp = Path(tempfile.mkdtemp(prefix='p3scale_'))
binary = P.build_binary()
w = (1 << np.arange(st.B.O)).astype(np.int64)
def masks_of(Y): return ((Y > 0).astype(np.int64) @ w)
def count(prob, extra):
    r = subprocess.run([binary, str(prob), '-t', '2', '-countlists'] + extra, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return {l.split()[1]: (int(l.split()[2]), int(l.split()[3])) for l in r.stdout.splitlines() if l.startswith('COUNT')}
ok = True
# ---- 1. filters
nB, nD = 600000, 300000
YB = st.B.sample(rows[1], nB, rng); YD = st.B.sample(rows[2], nD, rng)
mB, mD = masks_of(YB), masks_of(YD)
mB.astype(np.uint64).tofile(tmp / 'b.bin'); mD.astype(np.uint64).tofile(tmp / 'd.bin')
FB, FD = st.B.values(YB), st.B.values(YD)
psdB = (st.psd(FB) <= 2 * v + 0.5).all(1); psdD = (st.psd(FD) <= 4 * v + 0.5).all(1)
def canonB(m): return P.is_canonical(st, 0, int(m), 0) if False else None
def pm(perm, m): return sum(1 << i for i in range(st.B.O) if (m >> int(perm[i])) & 1)
permsB = [st.perm_of(st.B, t) for t in st.mult_reps[1:]]
revp = st.perm_of(st.B, st.neg)
def cB(m): return all(pm(p, m) >= m for p in permsB)
def cD(m): return pm(revp, m) >= m
t0 = time.time()
canB = np.array([cB(int(m)) for m in mB[:200000]]); canB_all = None
# canonical test on a 200k subset for the canonical-only comparison (python is slow); PSD+canonical on the PSD survivors of all
survB = [int(m) for m in mB[psdB]]
canB_surv = [cB(m) for m in survB]
canD = np.array([cD(int(m)) for m in mD])
survD_can = int(np.sum(psdD & canD))
print(f'python side {time.time()-t0:.1f}s', flush=True)
for label, canon, psd in (('PSD only', False, True), ('canonical only (B: first 200k masks)', True, False), ('PSD + canonical', True, True)):
    P.write_problem(st, rows, tmp / 'p.txt', canonB=canon, canonD=canon, psd=(False, psd, psd))
    if label.startswith('canonical only'):
        mB[:200000].astype(np.uint64).tofile(tmp / 'b200.bin')
        c = count(tmp / 'p.txt', ['-blist', str(tmp / 'b200.bin'), '-dlist', str(tmp / 'd.bin'), '-alist', str(tmp / 'b200.bin')]) if False else None
    cb = count(tmp / 'p.txt', ['-blist', str(tmp / ('b200.bin' if label.startswith('canonical only') else 'b.bin')), '-dlist', str(tmp / 'd.bin')])
    if label == 'PSD only':
        wantB, wantD = int(psdB.sum()), int(psdD.sum())
    elif label.startswith('canonical only'):
        wantB, wantD = int(canB.sum()), int(canD.sum())
    else:
        wantB, wantD = int(sum(canB_surv)), survD_can
    gotB, gotD = cb['B'][1], cb['D'][1]
    good = (gotB, gotD) == (wantB, wantD)
    ok &= good
    print(f'  filters [{label}]: engine kept B {gotB} D {gotD}; numpy B {wantB} D {wantD}: {"ok" if good else "MISMATCH"}', flush=True)
# ---- 2. DFS leaf counts on a slice (static task assignment: with -t 8 and slice s/64 only one thread works)
def tasks(blk, r, want):
    reach = None
    O = blk.O; sizes = blk.sizes
    tail = np.concatenate([np.cumsum(sizes[::-1])[::-1], [0]])
    W = blk._dp() if blk._W is None else blk._W
    def feasible(k, rem): return abs(rem) <= tail[k] and 0 <= rem + v <= 2 * v and W[k, rem + v] > 0
    for d in range(1, O - 1):
        out = []
        def rec(k, rem):
            if k == d: out.append((rem)); return
            for s in (1, -1):
                r2 = rem - s * int(sizes[k])
                if feasible(k + 1, r2): rec(k + 1, r2)
        rec(0, r)
        if len(out) >= want: break
    return d, out, W
P.write_problem(st, rows, tmp / 'p.txt', canonB=False, canonD=False, psd=(False, False, False))
for kind, blk, r, key in (('B', st.B, rows[1], 'B'), ('D', st.B, rows[2], 'D')):
    d, rems, W = tasks(blk, r, 64 * 8)
    sl = 5
    want = int(sum(W[d, rem + v] for i, rem in enumerate(rems) if i % 64 == sl))
    c = subprocess.run([binary, str(tmp / 'p.txt'), '-t', '8', '-countlists', '-slice', f'{sl}/64'], capture_output=True, text=True)
    got = int([l for l in c.stdout.splitlines() if l.startswith(f'COUNT {key} ')][0].split()[2])
    good = got == want
    ok &= good
    print(f'  DFS leaves {kind}, slice {sl}/64 (prefix depth {d}, {len(rems)} tasks): engine {got}, dynamic program {want}: {"ok" if good else "MISMATCH"}', flush=True)
print('ALL OK' if ok else 'FAILED'); sys.exit(0 if ok else 1)

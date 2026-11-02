"""Planted-solution test at the full dimensions of v = 239 (17 classes, 18 + 35 orbits, rows (15, 1, 27)).

A random triple (A0, B0, D0) is drawn from the true pools; the target T = F_A0 + 2 F_B0 + F_D0 is imposed (T != 0, so
the PSD and canonical filters are off: they assume a zero target).  A is enumerated in full (24310 members); the B and D
lists are random subsets (override files) containing the planted members.  The join must return exactly the planted
triple, for every split-class count ns = 1, 2, 3 and for 1 and 3 shards.   python3.11 p3plant.py [--nb N] [--nd N]"""
import argparse, sys, tempfile, time
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import p3run as P
from p3struct import Structure

ap = argparse.ArgumentParser()
ap.add_argument('--nb', type=int, default=40000); ap.add_argument('--nd', type=int, default=600000)
ap.add_argument('--seed', type=int, default=3); ap.add_argument('--threads', type=int, default=4)
a = ap.parse_args()
st = Structure(239, [10]); rows = (15, 1, 27)
rng = np.random.default_rng(a.seed)
w = (1 << np.arange(st.B.O)).astype(np.int64)
def masks(Y): return ((Y > 0).astype(np.int64) @ w)
YA0 = st.A.sample(rows[0], 1, rng); YB0 = st.B.sample(rows[1], 1, rng); YD0 = st.B.sample(rows[2], 1, rng)
tgt = st.A.values(YA0)[0] + 2 * st.B.values(YB0)[0] + st.B.values(YD0)[0]      # the engine solves F_A + 2 F_B + F_D = T
print('planted target', tgt.tolist())
mA0 = int(((YA0 > 0).astype(np.int64) @ (1 << np.arange(st.A.O)))[0]); mB0 = int(masks(YB0)[0]); mD0 = int(masks(YD0)[0])
tmp = Path(tempfile.mkdtemp(prefix='p3plant_'))
mb = np.concatenate([masks(YB0), masks(st.B.sample(rows[1], a.nb, rng))]).astype(np.uint64)
md = np.concatenate([masks(YD0), masks(st.B.sample(rows[2], a.nd, rng))]).astype(np.uint64)
mb.tofile(tmp / 'b.bin'); md.tofile(tmp / 'd.bin')
P.write_problem(st, rows, tmp / 'prob.txt', canonB=False, canonD=False, psd=(False, False, False), seed=5, tgt=tgt)
ok = True
for ns, shards in ((1, 1), (2, 1), (3, 1), (2, 3)):
    got = set(); t = time.time(); pairs = 0
    for sh in range(shards):
        r = P.run_p3(tmp / 'prob.txt', threads=a.threads, ns=ns, shard=(sh, shards) if shards > 1 else None, blist=tmp / 'b.bin', dlist=tmp / 'd.bin')
        got |= set(r['hits']); pairs += int(r['done']['pairs'])
        assert int(r['done']['badmatches']) == 0
    good = got == {(mA0, mB0, mD0)}
    ok &= good
    print(f'ns={ns} shards={shards}: lists {r["lists"]} pairs {pairs:.3e} hits {len(got)} planted found and alone: {good} [{time.time()-t:.1f}s]', flush=True)
# negative control: the same lists with a perturbed target must give no hit
tgt2 = tgt.copy(); tgt2[3] += 4
P.write_problem(st, rows, tmp / 'prob2.txt', canonB=False, canonD=False, psd=(False, False, False), seed=5, tgt=tgt2)
r = P.run_p3(tmp / 'prob2.txt', threads=a.threads, ns=2, blist=tmp / 'b.bin', dlist=tmp / 'd.bin')
neg = len(r['hits']) == 0
ok &= neg
print(f'negative control (target + 4 in class 3): hits {len(r["hits"])} (must be 0): {neg}')
print('ALL OK' if ok else 'FAILED'); sys.exit(0 if ok else 1)

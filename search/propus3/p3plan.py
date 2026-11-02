"""Cost preview of the exact reduced search for a structure (v, H): pool sizes, sampled PSD retention, expected list sizes
after the lossless filters, pair count, enumeration leaves, memory, and core-seconds at the rates measured on the Apple M4
(DFS 42-43 ns/leaf, two passes; join 3.4-10 ns/pair; D entries 16 bytes).  Retention is a uniform-sample estimate (validated
against exact engine counts on 532 and 876 in psd_retention_check.py).
   python3.11 p3plan.py V GEN [--samples 100000] [--rows a,b,d]"""
import argparse, sys, math
import numpy as np
sys.path.insert(0, '.')
from p3struct import Structure
import p3run as P

ap = argparse.ArgumentParser()
ap.add_argument('v', type=int); ap.add_argument('gen', type=int)
ap.add_argument('--samples', type=int, default=100000); ap.add_argument('--rows', default=None)
a = ap.parse_args()
st = Structure(a.v, [a.gen]); v = a.v
rng = np.random.default_rng(1)
print(f'v={v} H=<{a.gen}> |H|={len(st.H)} OA={st.A.O} OB={st.B.O} C={st.C} |Mult|={len(st.mult_reps)} |Gamma|={P.gamma_order(st)}')
pats = st.feasible_patterns()
if a.rows:
    want = tuple(int(x) for x in a.rows.split(',')); pats = [p for p in pats if p['rows'] == want]
tot = dict(leaves=0, pairs=0, mem=0, cs=0)
for p in pats:
    ra, rb, rd = p['rows']; pa, pb, pd = p['pools']
    rho = []
    for blk, r, b in ((st.A, ra, 4 * v), (st.B, rb, 2 * v), (st.B, rd, 4 * v)):
        Y = blk.sample(r, a.samples, rng)
        rho.append(float((st.psd(blk.values(Y)) <= b + 0.5).all(1).mean()))
    Ap = pa * rho[0]; Bp = pb * rho[1] / len(st.mult_reps); Dp = pd * rho[2] / 2
    pairs = Ap * Bp
    leaves = 2 * (pb + pd)
    cs = leaves * 43e-9 + pairs * 6e-9 + Dp * 15e-9
    print(f'  rows {p["rows"]}: pools A {pa:.3e} B {pb:.3e} D {pd:.3e}; PSD retention A {rho[0]:.4f} B {rho[1]:.5f} D {rho[2]:.4f}\n'
          f'     lists A\' {Ap:.3e} B\' {Bp:.3e} D\' {Dp:.3e}; pairs {pairs:.3e}; DFS leaves (2 passes) {leaves:.3e}; D memory {16 * Dp / 1e9:.1f} GB; '
          f'M4 core-seconds ~{cs:.0f} ({cs / 60:.1f} core-min)')
    tot['leaves'] += leaves; tot['pairs'] += pairs; tot['mem'] = max(tot['mem'], 16 * Dp / 1e9); tot['cs'] += cs
print(f'  TOTAL over patterns: leaves {tot["leaves"]:.3e}, pairs {tot["pairs"]:.3e}, peak D memory {tot["mem"]:.1f} GB, M4 core-seconds ~{tot["cs"]:.0f} ({tot["cs"] / 3600:.2f} core-h)')

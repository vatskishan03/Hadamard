"""Exhaustive reduced search of one structure (v, H) over all feasible row patterns.
  python3.11 p3sweep.py V GEN [--threads T] [--ns K] [--shards N] [--hits FILE]
Prints, per row pattern, the pool sizes, the filtered list sizes, the pairs probed, and the number of reduced hits
(one per Gamma-orbit for generic orbits; the group order is printed).  Every hit is re-checked by direct PAF."""
import argparse, sys, tempfile, time
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import p3run as P
from p3struct import Structure

ap = argparse.ArgumentParser()
ap.add_argument('v', type=int); ap.add_argument('gen', type=int)
ap.add_argument('--threads', type=int, default=4); ap.add_argument('--ns', type=int, default=2)
ap.add_argument('--shards', type=int, default=1); ap.add_argument('--hits', default=None)
ap.add_argument('--maxpool', type=float, default=1e12)
a = ap.parse_args()
st = Structure(a.v, [a.gen])
print(f'v={a.v} H=<{a.gen}> order {len(st.H)} OA={st.A.O} OB={st.B.O} C={st.C} |Mult|={len(st.mult_reps)} |Gamma|={P.gamma_order(st)}', flush=True)
tmp = Path(tempfile.mkdtemp(prefix='p3sweep_'))
tot = 0
for p in st.feasible_patterns():
    rows = p['rows']
    if max(p['pools']) > a.maxpool:
        print(f'  rows {rows}: pools {p["pools"]} skipped (> maxpool)'); continue
    pf = tmp / 'prob.txt'
    P.write_problem(st, rows, pf, canonB=True, canonD=True, psd=(True, True, True), seed=a.v)
    hits = set(); t = time.time(); lists = None; pairs = 0
    for sh in range(a.shards):
        r = P.run_p3(pf, threads=a.threads, ns=a.ns, shard=(sh, a.shards) if a.shards > 1 else None, hits=a.hits)
        hits |= set(r['hits']); lists = r['lists']; pairs += int(r['done']['pairs'])
    for h in hits:
        assert P.exact_check(st, *P.masks_to_seqs(st, *h))
    tot += len(hits)
    print(f'  rows {rows}: pools {p["pools"]}  lists(A,B,D) {lists}  pairs {pairs:.3e}  reduced hits {len(hits)}  [{time.time()-t:.1f}s]', flush=True)
print(f'total reduced hits (Gamma-orbit representatives, up to the translation symmetry of Fix(H)): {tot}')

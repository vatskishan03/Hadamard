"""Write the problem file read by p3 for one structure and row pattern.
   python3.11 p3problem.py V GEN --rows a,b,d --out PROBLEM [--noreduce] [--nopsd] [--seed N]
Default: canonical B (under Z_v^*/H) and canonical D (under reversal of D), PSD bounds on A, B, D (all lossless)."""
import argparse, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import p3run as P
from p3struct import Structure

ap = argparse.ArgumentParser()
ap.add_argument('v', type=int); ap.add_argument('gen', type=int)
ap.add_argument('--rows', required=True); ap.add_argument('--out', required=True)
ap.add_argument('--noreduce', action='store_true'); ap.add_argument('--nopsd', action='store_true')
ap.add_argument('--seed', type=int, default=None)
a = ap.parse_args()
st = Structure(a.v, [a.gen])
rows = tuple(int(x) for x in a.rows.split(','))
feas = {p['rows']: p['pools'] for p in st.feasible_patterns()}
if rows not in feas:
    sys.exit(f'row pattern {rows} is not feasible for this structure; feasible: {sorted(feas)}')
P.write_problem(st, rows, a.out, canonB=not a.noreduce, canonD=not a.noreduce, psd=(not a.nopsd,) * 3, seed=a.seed or a.v)
print(f'wrote {a.out}: v={a.v} |H|={len(st.H)} OA={st.A.O} OB={st.B.O} C={st.C} rows={rows} pools={feas[rows]} |Mult|={len(st.mult_reps)} |Gamma|={P.gamma_order(st)}')

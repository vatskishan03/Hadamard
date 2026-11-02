"""Write p3 problem files (same code path as propus239/p3problem.py: p3run.write_problem, canonical B and D, PSD bounds on all
three blocks, seed = v) for one structure (v, H = <gens>, any number of generators) and any number of row patterns.

  python3 -B mkprob.py V GENS --rows a,b,d [--rows ...] --outdir DIR --tag TAG
  GENS = comma separated generators of H in Z_v^* (p3problem.py itself only takes a single generator).
Files: DIR/TAG_r<a>_<b>_<d>.txt .  The tools are imported read-only from propus239 (run with python3 -B: no .pyc is written there)."""
import argparse, sys
from pathlib import Path
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))     # search/propus3
import p3run as P
from p3struct import Structure

ap = argparse.ArgumentParser()
ap.add_argument('v', type=int); ap.add_argument('gens')
ap.add_argument('--rows', action='append', required=True)
ap.add_argument('--outdir', required=True); ap.add_argument('--tag', required=True)
a = ap.parse_args()
st = Structure(a.v, [int(x) for x in a.gens.split(',')])
feas = {p['rows'] for p in st.feasible_patterns()}
Path(a.outdir).mkdir(parents=True, exist_ok=True)
for rs in a.rows:
    rows = tuple(int(x) for x in rs.split(','))
    if rows not in feas:
        sys.exit(f'row pattern {rows} is not feasible for v={a.v} gens={a.gens}')
    out = Path(a.outdir) / f'{a.tag}_r{rows[0]}_{rows[1]}_{rows[2]}.txt'
    P.write_problem(st, rows, out, canonB=True, canonD=True, psd=(True, True, True), seed=a.v)
    print(f'wrote {out}: v={a.v} |H|={len(st.H)} OA={st.A.O} OB={st.B.O} C={st.C} rows={rows} |Mult|={len(st.mult_reps)} |Gamma|={P.gamma_order(st)}')

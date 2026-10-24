"""Regression test of the correction of rowfast2.py (search/NOTES.md section 91): the vectorised row-sum count must not depend on which shapes the process built before.
Builds the shapes 575 (K = Z25 and Z5^2, H_K = 1), 1081 (Z23 x F47), 1139 (Z17 x F67) and 1207 (Z17 x F71), several times in alternation, in ONE process (each freed before the next), takes the slots of their top windows
(at most 25 each, the reference is slow for larger shapes) and compares `rowfast2.row_count_fast(sh, a, ms) > 0` with the reference `sk_pat.row_patterns_skew`.
The round-38 module kept its caches in module-level dictionaries keyed by id(sh); id() is reused after a shape is freed, so it used the block classes of another shape and reported feasible slots as infeasible.
usage (in search/skewplug): python3 scan38/test_rowfast2.py [path of another rowfast2.py to test, e.g. the round-38 original]   -> prints the number of mismatches (0 for the corrected module)"""
import sys, os, json, gc, random, importlib.util, warnings
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..'))
from pshape import shape_from_row
from sk_pat import m_splits, row_patterns_skew
from awin import alloc_window
if len(sys.argv) > 1:
    spec = importlib.util.spec_from_file_location('rowfast2_alt', sys.argv[1]); rf = importlib.util.module_from_spec(spec); spec.loader.exec_module(rf)
else:
    import rowfast2 as rf
rows = json.load(open(os.path.join(HERE, 'data', 'shapes38.json')))
want = [(575, 23, [25], 1), (575, 23, [5, 5], 1), (1081, 47, [23], 1), (1139, 67, [17], 1), (575, 23, [25], 1), (1207, 71, [17], 1), (575, 23, [5, 5], 1), (1081, 47, [23], 1), (575, 23, [25], 1)]
rng = random.Random(1)
bad = tot = 0
for (v, q, km, order) in want:
    r = next(x for x in rows if (x['v'], x['q'], x['kmods'], x['order']) == (v, q, km, order))
    sh = shape_from_row(r)
    win, emin, emax = alloc_window(sh, 5)
    slots = [(a, ms) for (a, em) in win for ms in m_splits(a, em)]
    rng.shuffle(slots)
    nb = 0
    for (a, ms) in slots[:25]:
        ref = bool(row_patterns_skew(sh, a, ms)); got = rf.row_count_fast(sh, a, ms) > 0
        tot += 1
        if ref != got: bad += 1; nb += 1
    print('shape', v, q, km, 'order', order, 'slots tested', min(25, len(slots)), 'mismatches', nb, flush=True)
    del sh; gc.collect()
print('TOTAL slots', tot, 'mismatches', bad)

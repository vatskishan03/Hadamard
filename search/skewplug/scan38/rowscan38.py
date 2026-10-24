"""EXACT all-slot row-sum scan (every allocation of every energy E, every m-split): is there ANY slot with an exact row-sum pattern (row_patterns_skew, a necessary condition)?
usage (in search/skewplug): python3 scan38/rowscan38.py SHAPES.json IDXLIST SECS  -> appends to runs/rowscan38.jsonl   status ALL_INFEASIBLE (EXACT: shape dead) | FEASIBLE (slot given) | PARTIAL (time cap)"""
import sys, os, json, time, warnings
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..'))      # scan38/ and the modules of search/skewplug
DATA = os.path.join(HERE, 'data')
from pshape import shape_from_row
from sk_pat import m_splits, row_patterns_skew
from rowfast2 import row_count_fast
from awin import alloc_window
os.makedirs('runs', exist_ok=True); rows = json.load(open(sys.argv[1])); idxs = [int(x) for x in sys.argv[2].split(',')]; secs = float(sys.argv[3])
for i in idxs:
    r = rows[i]; t0 = time.time()
    sh = shape_from_row(r)
    win, emin, emax = alloc_window(sh, 0, want_all=True)
    res = dict(idx=i, v=r['v'], q=r['q'], kmods=r['kmods'], order=r['order'], nu=r['nu'], nalloc=len(win), Emin=emin, Emax=emax)
    n = 0; found = None; status = 'ALL_INFEASIBLE'
    win.sort(key=lambda x: -sum(x[0]))
    for (a, em) in win:
        for ms in m_splits(a, em):
            if time.time() - t0 > secs: status = 'PARTIAL'; break
            n += 1
            try: p = row_count_fast(sh, a, ms)
            except Exception: p = None
            if p: found = (list(a), list(ms)); status = 'FEASIBLE'; break
        if status != 'ALL_INFEASIBLE': break
    res.update(status=status, slots_checked=n, found=found, secs=round(time.time() - t0, 1))
    open('runs/rowscan38.jsonl', 'a').write(json.dumps(res) + '\n')
    print(json.dumps(res), flush=True)

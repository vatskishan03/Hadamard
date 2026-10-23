"""validation battery: G-level brute force vs K-level brute force vs the two-level pipeline, ordered counts of skew GS quadruples (resume-safe -> val.jsonl)"""
import sys, os, json, time, warnings
warnings.filterwarnings('ignore')
from sk import SShape
from sk_exh import exhaust
from bf import gcount, kcount
SHAPES = []
for q in (3, 7, 11, 19, 23, 31, 43, 47): SHAPES.append(([3], q, []))
for q in (3, 7, 11, 19, 23, 31): SHAPES.append(([5], q, []))
for q in (3, 7, 11, 19, 23): SHAPES.append(([7], q, [[2]]))
for q in (3, 7, 11, 19): SHAPES.append(([9], q, [[4]]))
for q in (3, 7, 11, 19): SHAPES.append(([11], q, [[3]]))
for q in (3, 7, 11, 19): SHAPES.append(([13], q, [[3]]))
done = set()
if os.path.exists('val.jsonl'):
    for l in open('val.jsonl'):
        d = json.loads(l); done.add((tuple(d['K']), d['q'], json.dumps(d['hu'])))
for (K, q, hu) in SHAPES:
    key = (tuple(K), q, json.dumps(hu))
    if key in done: continue
    t0 = time.time(); sh = SShape(K, q, hu)
    rec = dict(K=K, q=q, hu=hu, G=sh.G)
    try: rec['gcount'] = gcount(K, q, hu)[0]
    except Exception as ex: rec['gcount'] = None; rec['gerr'] = str(ex)[:80]
    try: rec['kcount'] = kcount(sh)[0]
    except Exception as ex: rec['kcount'] = None; rec['kerr'] = str(ex)[:80]
    big = max([v for v in (rec['gcount'], rec['kcount']) if v is not None] or [0])
    if big > 3e8:
        rec['pipe'] = None; rec['status'] = 'G_vs_K_ONLY(count too large for the pipeline)' if rec['gcount'] == rec['kcount'] else 'MISMATCH'; rec['secs'] = round(time.time() - t0, 1)
        open('val.jsonl', 'a').write(json.dumps(rec) + '\n'); print(rec, flush=True); continue
    tot, sols = exhaust(sh, tmax=150)
    rec['pipe'] = tot['E_G']; rec['pipe_exhaustive'] = tot['exhaustive']; rec['allocs'] = tot['allocs']; rec['n1_G'] = tot['N1_G']
    if sols: rec['sol_check_G'] = bool(sh.g_target_skew_ok(sh.glevel_blocks(sols[0])))
    vals = [v for v in (rec['gcount'], rec['kcount'], rec['pipe']) if v is not None]
    rec['status'] = 'MATCH' if len(set(vals)) == 1 and rec['pipe_exhaustive'] else ('MISMATCH' if len(set(vals)) > 1 else 'PIPE_NOT_EXHAUSTIVE')
    rec['secs'] = round(time.time() - t0, 1)
    open('val.jsonl', 'a').write(json.dumps(rec) + '\n')
    print(rec, flush=True)

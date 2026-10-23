"""Aggregates production logs (runs/*.jsonl of prod_ks.py, prod_dedup.py, sk_worker.py): per worker the number of rounds, slots with level-1 instances, level-1 instances, instances given to level 2,
instances with a completion, seconds in level 1 and in level 2, level-2 instances per second, and the order of a hit.  For logs of prod_dedup.py also the number of NEW (distinct) classes and the duplicate factor.
usage (from search/skewplug):  python3 estimates/summarize_runs.py [GLOB ...]      default: runs/*.jsonl"""
import sys, json, glob, os
rows_out = []
for f in sorted(g for pat in (sys.argv[1:] or ['runs/*.jsonl']) for g in glob.glob(pat)):
    name = os.path.basename(f)[:-6]
    recs = []
    for l in open(f):
        try: recs.append(json.loads(l))
        except Exception: pass
    recs = [r for r in recs if 'al' in r]
    if not recs: continue
    new = sum(r.get('n1_new', 0) for r in recs)
    n1 = sum(r.get('n1', 0) for r in recs); l2 = sum(r.get('n_l2', 0) for r in recs); nok = sum(r.get('nok', 0) for r in recs)
    t1 = sum(r.get('t_l1', 0) for r in recs); t2 = sum(r.get('t_l2', 0) for r in recs)
    hit = [r['HIT'] for r in recs if 'HIT' in r]
    slots = len({(tuple(r['al']), tuple(r['ms'])) for r in recs if r.get('n1', 0) > 0})
    rows_out.append((name, len(recs), slots, n1, l2, nok, round(t1), round(t2), round(l2 / max(t1 + t2, 1e-9)), new, round(sum(r['n1'] for r in recs if 'n1_new' in r) / new, 1) if new else 0, 'HIT ' + str(hit[0]['order']) if hit else ''))
print(f'{"worker":16} {"rounds":>6} {"slots":>5} {"level1":>9} {"level2":>9} {"complet.":>8} {"t1 s":>6} {"t2 s":>6} {"l2/s":>6} {"new":>8} {"dup":>6}')
for r in rows_out: print(f'{r[0]:16} {r[1]:6d} {r[2]:5d} {r[3]:9d} {r[4]:9d} {r[5]:8d} {r[6]:6d} {r[7]:6d} {r[8]:6d} {r[9]:8d} {r[10]:6.1f} {r[11]}')

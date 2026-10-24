"""exh_table.py -- the table of the exhaustive plug-in verdicts of search/NOTES.md section 92 from the recorded logs of exh_sweep.py (../data/gs/exh_small.jsonl, exh_r3.jsonl, exh_alive.jsonl; or the files named).
A shape is EXHAUSTIVE-empty when gs_exh.exhaust ran over all allocations with no skipped allocation and no flagged level-2 instance ('exhaustive' true) and found E_G = 0 ordered quadruples.
Prints one line per distinct shape (mode, n, r, K, H) and the totals per p (n = 2p); shapes that did not finish are listed with the reason (UNKNOWN).
usage: exh_table.py [FILE.jsonl ...]"""
import sys, os, json, collections
HERE = os.path.dirname(os.path.abspath(__file__))
files = sys.argv[1:] or [os.path.join(HERE, '..', 'data', 'gs', f) for f in ('exh_small.jsonl', 'exh_r3.jsonl', 'exh_alive.jsonl')]
best = {}
for f in files:
    for l in open(f):
        d = json.loads(l)
        k = d['skey']
        ok = d.get('exhaustive') is True and d.get('E_G') == 0
        if k not in best or (ok and not best[k][0]): best[k] = (ok, d)
done = collections.defaultdict(list); undone = collections.defaultdict(list)
for k, (ok, d) in best.items():
    p = d['n'] // 2
    row = (d['orb'], d['K'], d['r'], d['hg'], d.get('secs'))
    if ok: done[p].append(row)
    else: undone[p].append(row + ('error ' + d['error'][:40] if d.get('error') else 'not exhaustive (time limit)',))
for p in sorted(set(done) | set(undone)):
    print('p = %d (n = 2p = %d): %d shapes EXHAUSTIVE-empty, %d not finished' % (p, 2 * p, len(done[p]), len(undone[p])))
    for row in sorted(done[p], key=lambda r: (r[0], str(r[1]), r[2], str(r[3]))): print('   EMPTY   %2d orbits  K = %s  r = %d  H = %s   (%.1f s)' % (row[0], row[1], row[2], row[3], row[4]))
    for row in sorted(undone[p], key=lambda r: (r[0], str(r[1]), r[2], str(r[3]))): print('   UNKNOWN %2d orbits  K = %s  r = %d  H = %s   (%s)' % (row[0], row[1], row[2], row[3], row[5]))
print('total: %d distinct shapes EXHAUSTIVE-empty, %d not finished' % (sum(len(v) for v in done.values()), sum(len(v) for v in undone.values())))

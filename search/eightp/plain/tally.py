import glob, ast, sys, os
# usage: tally.py [DIR]   DIR holds spin/ and plain/ with the logs (default ../data = the recorded logs; use runs/-directories after a re-run)
D = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
def rows(f):
    for l in open(f):
        l = l.strip()
        if l.startswith('{'):
            try: yield ast.literal_eval(l)
            except Exception: pass
print('== spin d=2/d=4 (spinD)')
for f in sorted(glob.glob(D + '/spin/D*.log')):
    hdr = open(f).readline().strip(); rs = list(rows(f))
    live = [r for r in rs if r.get('L') != 0]
    print(f.split('/')[-1], '|', hdr[:60], '| done pairs', len(rs), 'live', len(live), 'maxOH_live', max([r['OH'] for r in live], default=0), 'count', sum(r['count'] for r in rs), 'secs', round(sum(r['secs'] for r in rs)))
print('== spin d=3')
for f in sorted(glob.glob(D + '/spin/v*.log')):
    hdr = open(f).readline().strip(); rs = list(rows(f))
    live = [r for r in rs if r.get('NA', 0) > 0]
    print(f.split('/')[-1], '|', hdr[:60], '| done pairs', len(rs), 'live', len(live), 'maxOH_live', max([r['OH'] for r in live], default=0), 'count', sum(r['count'] for r in rs), 'secs', round(sum(r['secs'] for r in rs)))
print('== plain (corrected solver r3/r4)')
for f in sorted(glob.glob(D + '/plain/r[34]_*.log')):
    rs = list(rows(f)); st = {}
    for r in rs: st[r['status']] = st.get(r['status'], 0) + 1
    err = any('Error' in l for l in open(f))
    print(f.split('/')[-1], 'shapes done', len(rs), st, 'count', sum(r.get('count_rho_plus', 0) or 0 for r in rs), 'ERROR' if err else '')

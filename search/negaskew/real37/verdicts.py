import json, glob, re
rows = {}
for f in sorted(glob.glob('runs/q*.jsonl')):
    for l in open(f):
        try: r = json.loads(l)
        except Exception: continue
        key = r['name']
        # keep the best status (EMPTY/HIT beats TIMEOUT)
        if key in rows and rows[key]['status'] in ('EMPTY', 'HIT') and r['status'] not in ('EMPTY', 'HIT'): continue
        rows[key] = dict(r, file=f)
out = []
for k, r in sorted(rows.items(), key=lambda kv: (kv[1]['v'], kv[0])):
    st = r.get('stats') or ''
    m = re.search(r'l1_instances=(\d+)', st); l1 = m.group(1) if m else '-'
    m = re.search(r'l2_instances=(\d+)', st); l2 = m.group(1) if m else '-'
    m = re.search(r'probes=(\d+)', st); pr = m.group(1) if m else '-'
    out.append('%-6s v=%-5s %-34s q=%-3s sec=%-8s l1=%-10s l2=%-9s probes=%-12s %s' % (r['status'], r['v'], k, r['q'], r['sec'], l1, l2, pr, (r.get('lastprog') or '')[:60] if r['status'] != 'EMPTY' else ''))
print('\n'.join(out)); print(len(out), 'jobs;', sum(1 for r in rows.values() if r['status'] == 'EMPTY'), 'EMPTY;', sum(1 for r in rows.values() if r['status'] == 'HIT'), 'HIT')

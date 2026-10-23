"""merge est/NAME.ext.jsonl (feasible slots with full logE) into est/NAME.jsonl (keeps the original lines; adds new (al, ms)); usage: est_merge.py NAME [MINLOGE]"""
import sys, json, os, shutil
name = sys.argv[1]; mn = float(sys.argv[2]) if len(sys.argv) > 2 else -99
p = f'est/{name}.jsonl'; e = f'est/{name}.ext.jsonl'
if not os.path.exists(p + '.orig'): shutil.copy(p, p + '.orig')
rows = [json.loads(l) for l in open(p + '.orig')]
have = {(tuple(r['al']), tuple(r['ms'])) for r in rows}
add = 0
import glob
for l in (l for f in sorted(glob.glob(f'est/{name}.ext*.jsonl')) for l in open(f)):
    r = json.loads(l)
    if 'logE' not in r or r['logE'] < mn: continue
    k = (tuple(r['al']), tuple(r['ms']))
    if k in have: continue
    have.add(k); rows.append(dict(al=r['al'], ms=r['ms'], logN1=r['logN1'], logE=r['logE'], logrho=r['logrho'])); add += 1
open(p + '.tmp', 'w').write(''.join(json.dumps(r) + '\n' for r in rows)); os.replace(p + '.tmp', p)
rows.sort(key=lambda r: -r['logE'])
print(name, 'orig', len(have) - add, 'added', add, 'top', [(r['al'], r['ms'], r['logE']) for r in rows[:6]])

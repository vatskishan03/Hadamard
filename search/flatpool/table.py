#!/usr/bin/env python3
"""table.py : predicted flat-pool cost for the measured structures (tg_<v>_<u>.json) -> table.md / table.json"""
import glob, json, math, sys
import predict
FR = {739: (9, 33.5, 22.7), 955: (19, 2.9, 28.8), 1103: (19, 21.1, 15.5), 1223: (13, 33.7, 32.4), 1257: (22, 12.7, 37.8), 1291: (15, 30.2, 29.3),
      1667: (17, 32.0, 37.5), 1713: (38, 5.4, 32.1), 1787: (19, 30.0, 36.1), 1893: (42, 4.4, 33.1), 2215: (34, 1.9, 50.1), 2971: (55, 14.0, 19.2)}
rows = []
for fn in sorted(glob.glob('tg_*_*.json')):
    if 'prelim' in fn: continue
    d = json.load(open(fn))
    if not d.get('rem'): continue
    row = dict(v=d['v'], ulen=d['ulen'], C=d['C'], pattern=d['pattern'])
    for key, mem, mode in (('3GB_chunk', 3.0, 'chunk'), ('3GB_part', 3.0, 'partition'), ('64GB_chunk', 64.0, 'chunk'), ('256GB', 256.0, 'chunk'), ('2TB', 2048.0, 'chunk')):
        dd, mod = predict.load(fn, mem)
        b = predict.scan(mod, mode=mode)
        if b is None: row[key] = None; continue
        tot, r = b
        row[key] = dict(total=tot, phi=r['phi'], f=r['f'], nneck=r['n'][0], Ecl=r['E_classes'], gen=r['gen'], join=r['join'][0], q=r['join'][1], s=r['join'][2])
    rows.append(row)
json.dump(rows, open('table.json', 'w'), indent=1)
def l(x): return "%.1f" % math.log10(x) if x else "-"
out = ["| n=4v | v | \\|U\\| | C | row sums | triage log10 E | triage log10 W | flat: phi, pool frac, necklaces/role, E_classes (2 TB) | log10 cpu-s: 3 GB part./chunk | 64 GB | 256 GB | 2 TB |", "|---|---|---|---|---|---|---|---|---|---|---|---|"]
for r in rows:
    fr = FR.get(r['v'], (None, None, None)); t = r['2TB']
    s2 = "phi %.2f, f %.2g, %.1e, E_cl %.1e" % (t['phi'], t['f'], t['nneck'], t['Ecl']) if t else '-'
    g = lambda k: l(r[k]['total']) if r[k] else '-'
    out.append("| %d | %d | %d | %d | %s | %.1f | %.1f | %s | %s / %s | %s | %s | %s |" % (4 * r['v'], r['v'], r['ulen'], r['C'], tuple(r['pattern']), fr[1], fr[2], s2, g('3GB_part'), g('3GB_chunk'), g('64GB_chunk'), g('256GB'), g('2TB')))
open('table.md', 'w').write("\n".join(out) + "\n")
print("\n".join(out))

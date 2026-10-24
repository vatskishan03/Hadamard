"""summarise est/*.ext.jsonl + hdr.json -> runs/summary38.json / runs/SUMMARY38.txt (ranked).  usage (in search/skewplug): python3 scan38/summ38.py [SHAPES.json]   (default scan38/data/shapes38.json)"""
import json, glob, os, sys, math
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data'); os.makedirs('runs', exist_ok=True)
shapes = json.load(open(sys.argv[1] if len(sys.argv) > 1 else DATA + '/shapes38.json'))
out = []
for hf in sorted(glob.glob('est/S*.hdr.json')):
    h = json.load(open(hf)); name = h['name']; idx = int(name[1:5])
    r = shapes[idx]
    rec = dict(window=h.get('window', 'top'), idx=idx, name=name, v=h['v'], q=h['q'], kmods=h['kmods'], nu=h['nu'], order=h['order'], ktype=r.get('ktype'), partial_enum=r.get('partial'), decided=r.get('decided'),
               status=h.get('status'), pool=h.get('pool'), processed=h.get('processed', 0), coverage=h.get('coverage'), nalloc=h.get('nalloc'), Emin=h.get('Emin'), Emax=h.get('Emax'), secs=h.get('secs'))
    ef = f'est/{name}.ext.jsonl'; feas = []; nopat = 0; nerr = 0; npr = 0
    if os.path.exists(ef):
        for l in open(ef):
            try: x = json.loads(l)
            except Exception: continue
            npr += 1
            if 'logE' in x: feas.append(x)
            elif 'nopat' in x: nopat += 1
            else: nerr += 1
    rec.update(processed=npr, nfeas=len(feas), nnopat=nopat, nerr=nerr, nunchecked=sum(1 for x in feas if x.get('unchecked')))
    pl = h.get('pool')
    if pl and pl == pl: rec['coverage'] = float('%.3g' % (npr / pl))
    elif h.get('window', 'top') == 'top' and r.get('pool5_ub'): rec['coverage'] = float('%.3g' % (npr / r['pool5_ub'])); rec['pool'] = r['pool5_ub']
    else: rec['coverage'] = None
    if feas:
        feas.sort(key=lambda x: -x['logE'])
        b = feas[0]
        rec.update(best_logE=b['logE'], best_logrho=b['logrho'], best_logN1=b['logN1'], best_slot=[b['al'], b['ms']])
        rec['top40'] = round(sum(10 ** x['logE'] for x in feas[:40]), 3)
        both = [x for x in feas if x['logE'] >= 0 and x['logrho'] >= -8.5]
        near = [x for x in feas if x['logE'] >= -1 and x['logrho'] >= -9.5]
        rec['n_both'] = len(both); rec['n_near'] = len(near)
        nuS = h['nu']; NN = r['m']
        def tbits(x): return nuS * (1 - min(x['al'][1:]) / NN)
        for x in feas: x['tb'] = tbits(x)
        f16 = [x for x in feas if x['tb'] <= 16.5]
        if f16:
            m16 = max(f16, key=lambda x: min(x['logE'], x['logrho'] + 8.5)); rec['margin16'] = round(min(m16['logE'], m16['logrho'] + 8.5), 2); rec['margin16_slot'] = [m16['al'], m16['ms'], m16['logE'], m16['logrho'], round(m16['tb'], 1)]
        else: rec['margin16'] = None
        mg = max(feas, key=lambda x: min(x['logE'], x['logrho'] + 8.5))
        rec['tbits_at_margin_slot'] = round(mg['tb'], 1)
        rec['margin'] = round(min(mg['logE'], mg['logrho'] + 8.5), 2); rec['margin_slot'] = [mg['al'], mg['ms'], mg['logE'], mg['logrho']]
        rec['best_rho'] = max(x['logrho'] for x in feas)
        rec['best_rho_slot_logE'] = max(feas, key=lambda x: x['logrho'])['logE']
        # best logrho among slots with logE >= -1
        c = [x for x in feas if x['logE'] >= -1]
        rec['best_rho_given_E'] = max((x['logrho'] for x in c), default=None)
    out.append(rec)
json.dump(out, open('runs/summary38.json', 'w'))
def key(r): return -r.get('margin', -999)
out.sort(key=key)
with open('runs/SUMMARY38.txt', 'w') as f:
    f.write('%-5s %-6s %-4s %-14s %-4s %-5s %-7s %-6s %-7s %-7s %-7s %-7s %-6s %s\n' % ('idx', 'v', 'q', 'K', 'nu', 'ord', 'pool', 'cov', 'nfeas', 'bestE', 'rho@E', 'top40', 'margin', 'status/notes'))
    for r in out:
        f.write('%-5s %-6s %-4s %-14s %-4s %-5s %-7s %-6s %-7s %-7s %-7s %-7s %-6s %s %s\n' % (r['idx'], r['v'], r['q'], 'x'.join(map(str, r['kmods'])), r['nu'], r['order'], r['pool'], r.get('coverage'), r.get('nfeas'), r.get('best_logE'), r.get('best_logrho'), r.get('top40'), r.get('margin'), r['status'], ('DECIDED' if r['decided'] else '')))
print(len(out), 'shapes summarised')

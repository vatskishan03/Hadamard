"""final table (in search/skewplug): python3 scan38/table38.py  -> runs/TABLE38.txt, runs/TABLE38_decided.txt, runs/table38.json   (merges scan38/data/shapes38.json, runs/summary38.json, runs/rowscan38.jsonl, lemma/manual overrides)"""
import json, os, sys, math
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data'); os.makedirs('runs', exist_ok=True)
from math import isqrt
shapes = json.load(open(DATA + '/shapes38.json'))
summ = {}; summ_top = {}
for s_ in json.load(open('runs/summary38.json')):
    if s_.get('window', 'top') == 'lo': summ[s_['idx']] = s_
    else: summ_top[s_['idx']] = s_
for i_, s_ in summ_top.items():
    if i_ not in summ: summ[i_] = s_
rs = {}
if os.path.exists('runs/rowscan38.jsonl'):
    for l in open('runs/rowscan38.jsonl'):
        try: r = json.loads(l); rs[r['idx']] = r
        except Exception: pass
def two_squares(n):
    for a in range(0, isqrt(n) + 1):
        b2 = n - a * a; b = isqrt(b2)
        if b * b == b2: return True
    return False
# shapes decided in the prose of NOTES 88/89/90 or r34-r37 STATUS (matched by hand on (v, q, K, order, nu)); the automatic match only sees job names
PROSE = {(201, 3, (67,), 11, 7): 'r38: Z67 x F3 <9> EXHAUSTIVELY EMPTY (r37 ic 0-19 + r38 resume ic 19-53, identical lists; 2.28e9 level-1, 2.5e8 level-2 instances, 0 solutions)', (575, 23, (5, 5), 3, 9): 'r34 sec.C / r37: Z5^2 x F23 Singer-3 EXHAUSTIVELY EMPTY', (575, 23, (5, 5), 5, 9): 'r34 sec.C (r33): Z5^2 x F23 order 5 EXHAUSTIVELY EMPTY',
         (1183, 7, (13, 13), 21, 9): 'r37: Z13^2 x F7 Singer-21 EXHAUSTIVELY EMPTY (4905 level-1 instances)', (1083, 3, (19, 19), 45, 9): 'r37: Z19^2 x F3 Singer-45 EXHAUSTIVELY EMPTY (no A list)',
         (1859, 11, (13, 13), 21, 9): 'r34: Z13^2 x F11 Singer-21 EXACTLY EMPTY (weights)', (2375, 19, (5, 5, 5), 31, 5): 'r34: Z5^3 x F19 Singer-31 EXACTLY EMPTY (weights)',
         (1029, 3, (7, 7, 7), 171, 3): 'r34: Z7^3 x F3 Singer-171 EXHAUSTIVELY EMPTY (brute force)', (2401, 7, (7, 7, 7), 171, 3): 'r34: Z7^3 x F7 Singer-171 EXHAUSTIVELY EMPTY (brute force)'}
import math
def shortfall(c, p90=False):
    # best-slot shortfall of a uniform random sample of fraction c of the slots, from bootstrap on the complete controls 539 / 781 / 1243 (median, p90 over 400 resamples)
    if c is None: return None
    if c >= 1: return 0.0
    c = max(c, 1e-9)
    pts = [(1.0, 0.0, 0.0), (0.3, 0.12, 0.3), (0.1, 0.3, 0.95), (0.03, 0.8, 2.0), (0.01, 1.15, 2.4)]
    for (c1, m1, p1), (c2, m2, p2) in zip(pts, pts[1:]):
        if c >= c2:
            f = (math.log10(c1) - math.log10(c)) / (math.log10(c1) - math.log10(c2)); return (p1 + f * (p2 - p1)) if p90 else (m1 + f * (m2 - m1))
    extra = 0.5 * math.log10(0.01 / c)
    return (2.4 if p90 else 1.15) + extra
out = []
for r in shapes:
    i = r['idx']; s = summ.get(i, {})
    st = None; note = ''
    kt = tuple(r['kmods'])
    if r['decided']: st = 'DECIDED'; note = r['decided'][:70]
    if r.get('nalloc') == 0: st = 'EXACT_DEAD(identity)' if not st else st + '+identity-dead'
    p = PROSE.get((r['v'], r['q'], kt, r['order'], r['nu']))
    if p and not r['decided']: st = 'DECIDED(prose)'; note = p
    if not st:
        # lemmas of NOTES 15 / 88 and the r34 big-big filter (hk1_batch.jsonl: H_K = 1, every pattern dead)
        if r['order'] == 1 and (r['v'], r['q'], r['m']) in {(361, 19, 19), (529, 23, 23), (899, 31, 29), (961, 31, 31), (1269, 27, 47), (1763, 43, 41)}: st = 'EXACT_DEAD(r34 big-big filter, H_K = 1, NOTES 88)'
        elif r.get('Emax') == 1 and not two_squares(2 * r['m'] - 1): st = 'EXACT_DEAD(N_w=1 lemma, 2m-1 not a sum of two squares)'
        elif r['v'] in (1081, 1711) and r.get('Emax') == 2: st = 'EXACT_DEAD(N_w=2 lemma, NOTES 88)'
    if not st and i in rs and rs[i]['status'] == 'ALL_INFEASIBLE': st = 'EXACT_DEAD(no row-sum pattern in any slot)'
    if not st: st = 'SCANNED' if i in summ else 'NOT_SCANNED'
    rec = dict(idx=i, v=r['v'], q=r['q'], kmods=r['kmods'], ktype=r['ktype'], order=r['order'], nu=r['nu'], gens=r['gens'][:2], partial_enum=r['partial'], status=st, note=note,
               Emin=r.get('Emin'), Emax=r.get('Emax'), pool_ub=r.get('pool5_ub'))
    rec.update({k: s.get(k) for k in ('pool', 'processed', 'coverage', 'nfeas', 'best_logE', 'best_logrho', 'best_logN1', 'best_slot', 'top40', 'margin', 'margin_slot', 'margin16', 'margin16_slot', 'tbits_at_margin_slot', 'nunchecked', 'n_both', 'best_rho', 'scan_status') if k in s})
    if s: rec['scan_status'] = s.get('status'); rec['window'] = s.get('window', 'top')
    if rec.get('margin') is not None and rec.get('coverage') is None:
        rec['margin_adj_med'] = rec['margin']; rec['margin_adj_p90'] = rec['margin']      # lazily sampled windows (q = 3 low window, |K| > 150): coverage unknown, no adjustment
    if rec.get('margin') is not None and rec.get('coverage') is not None:
        rec['margin_adj_med'] = round(rec['margin'] + shortfall(rec['coverage']), 2); rec['margin_adj_p90'] = round(rec['margin'] + shortfall(rec['coverage'], True), 2)
    out.append(rec)
json.dump(out, open('runs/table38.json', 'w'))
sc = [o for o in out if o['status'] == 'SCANNED']
sc.sort(key=lambda o: (-(o.get('margin_adj_med') if o.get('margin_adj_med') is not None else -9999), -(o.get('best_logE') if o.get('best_logE') is not None else -9999)))
with open('runs/TABLE38.txt', 'w') as f:
    f.write('A. UNDECIDED shapes scanned (ranked by margin = max over slots of min(logE, logrho + 8.5); criterion: margin >~ 0).  controls: 473 +0.67/-7.0 (margin +0.67), 539 <(2,2)> +0.94/-7.87 (+0.63)\n')
    f.write('%-5s %-5s %-4s %-12s %-5s %-4s %-6s %-9s %-7s %-6s %-8s %-8s %-8s %-8s %-8s %-8s %-6s %-8s %s\n' % ('idx', 'v', 'q', 'K', 'ord', 'nu', 'E', 'pool', 'cover', 'nfeas', 'bestE', 'rho@E', 'margin', 'adjMed', 'adjP90', 'tbits', 'm16', 'top40', 'scan/notes'))
    for o in sc:
        f.write('%-5s %-5s %-4s %-12s %-5s %-4s %-6s %-9s %-7s %-6s %-8s %-8s %-8s %-8s %-8s %-6s %-8s %-8s %s%s\n' % (o['idx'], o['v'], o['q'], 'x'.join(map(str, o['kmods'])), o['order'], o['nu'], '%s-%s' % (o['Emin'], o['Emax']), ('%.3g' % o['pool']) if o.get('pool') else '-', ('%.3g' % o['coverage']) if o.get('coverage') is not None else '-', o.get('nfeas'), o.get('best_logE'), o.get('best_logrho'), o.get('margin'), o.get('margin_adj_med'), o.get('margin_adj_p90'), o.get('tbits_at_margin_slot'), o.get('margin16'), o.get('top40'), str(o.get('scan_status')) + '/' + str(o.get('window')), ' PARTIAL-ENUM' if o['partial_enum'] else ''))
    f.write('\nB. UNDECIDED shapes EXACTLY dead / not scanned\n')
    for o in out:
        if o['status'].startswith('EXACT_DEAD') and not o['status'].startswith('EXACT_DEAD(identity)') and not o['note'] or o['status'] == 'NOT_SCANNED':
            f.write('%s v=%s q=%s K=%s ord=%s nu=%s  %s\n' % (o['idx'], o['v'], o['q'], 'x'.join(map(str, o['kmods'])), o['order'], o['nu'], o['status']))
    f.write('\nC. counts: %s\n' % {s: sum(1 for o in out if o['status'] == s) for s in set(o['status'] for o in out)})
print(open('runs/TABLE38.txt').read()[:6000])

with open('runs/TABLE38_decided.txt', 'w') as f:
    f.write('DECIDED / EXACT-DEAD shapes of the enumeration (not scanned; verdicts of rounds 34-37 or exact filters of this round)\n')
    for o in sorted(out, key=lambda o: (o['v'], o['q'], o['order'])):
        if o['status'] != 'SCANNED' and o['status'] != 'NOT_SCANNED':
            f.write('%-5s v=%-5s q=%-3s K=%-10s ord=%-5s nu=%-3s %s  %s\n' % (o['idx'], o['v'], o['q'], 'x'.join(map(str, o['kmods'])), o['order'], o['nu'], o['status'], o['note']))

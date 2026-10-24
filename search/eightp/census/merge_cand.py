"""merge_cand.py -- the merged census files of ../data/census from the per-n outputs of census_a.py (search/NOTES.md section 92):
   cand8p.jsonl     = cens_<n>.jsonl for n = 1198, 1438, 1478, 1678, 2062, 2174 followed by censP_<n>.jsonl (cyclic and product multiplier groups, r >= 5)
   cand3only.jsonl  = the rows of cens3_<n>.jsonl (r >= 3) that are not in cens_<n>.jsonl (the shapes with r = 3);  cand3Ponly.jsonl likewise from cens3P and censP
   candC8p.jsonl    = censC_<n>.jsonl (mode c, bordered)
   est_all8p.jsonl  = est8p.jsonl + est3only.jsonl + est3Ponly.jsonl (the estimates of run_est.py on the live rows)
usage: merge_cand.py [DATADIR]      compares the recomputed files with the recorded ones and prints OK or the number of differing lines; merge_cand.py DATADIR --write OUTDIR writes them"""
import sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__))
D = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith('--') else os.path.join(HERE, '..', 'data', 'census')
NS = (1198, 1438, 1478, 1678, 2062, 2174)
def load(tag): return [l for n in NS for l in open(os.path.join(D, '%s_%d.jsonl' % (tag, n)))]
def key(l):
    d = json.loads(l); return json.dumps([d['mode'], d['n'], d['r'], d['K'], d['hg'], d['qspec']])
def read(f): return open(os.path.join(D, f)).readlines()
cens, censP, cens3, cens3P, censC = (load(t) for t in ('cens', 'censP', 'cens3', 'cens3P', 'censC'))
kc = {key(l) for l in cens}; kP = {key(l) for l in censP}
out = {'cand8p.jsonl': cens + censP, 'cand3only.jsonl': [l for l in cens3 if key(l) not in kc], 'cand3Ponly.jsonl': [l for l in cens3P if key(l) not in kP], 'candC8p.jsonl': censC,
       'est_all8p.jsonl': read('est8p.jsonl') + read('est3only.jsonl') + read('est3Ponly.jsonl')}
if '--write' in sys.argv:
    o = sys.argv[sys.argv.index('--write') + 1]; os.makedirs(o, exist_ok=True)
    for f, ls in out.items(): open(os.path.join(o, f), 'w').writelines(ls)
    print('written to', o); sys.exit()
for f, ls in out.items():
    rec = read(f)
    print('%-18s %4d rows  %s' % (f, len(ls), 'OK (identical to the recorded file)' if ls == rec else 'DIFFERS: %d recorded rows, %d recomputed, %d equal as sets' % (len(rec), len(ls), len(set(rec) & set(ls)))))
live = [l for l in out['cand8p.jsonl'] if json.loads(l)['status'] == 'live']
print('cand8p: %d rows, %d live' % (len(out['cand8p.jsonl']), len(live)))

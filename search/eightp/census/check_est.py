"""check_est.py -- re-estimates (run_est.py with the arguments of the recorded run: M = 700, nfull = 150, topk = 6) the live rows of ../data/census/cand8p.jsonl with n = N and compares the estimates with ../data/census/est8p.jsonl.
usage: check_est.py N [N ...]     (N = 1478 takes 10 s, 7 live shapes)       prints 'k of m identical'"""
import sys, os, json, subprocess
HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, '..', 'data', 'census'); R = os.path.join(HERE, 'runs'); os.makedirs(R, exist_ok=True)
ns = [int(x) for x in sys.argv[1:]]
rows = [l for l in open(os.path.join(D, 'cand8p.jsonl')) if json.loads(l)['n'] in ns]
inp = os.path.join(R, 'est_in.jsonl'); out = os.path.join(R, 'est_out.jsonl')
open(inp, 'w').writelines(rows)
if os.path.exists(out): os.remove(out)
subprocess.run([sys.executable, os.path.join(HERE, 'run_est.py'), inp, out, '700', '150', '6'], check=True, stdout=subprocess.DEVNULL)
a = {json.loads(l)['key']: json.loads(l)['est'] for l in open(out)}
b = {json.loads(l)['key']: json.loads(l)['est'] for l in open(os.path.join(D, 'est8p.jsonl')) if json.loads(l).get('est')}
same = sum(1 for k in a if b.get(k) == a[k])
print('%d of %d identical (logE of the best shapes: %s)' % (same, len(a), sorted(round(v['logE'], 2) for v in a.values())[-3:]))

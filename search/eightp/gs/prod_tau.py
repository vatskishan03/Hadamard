"""prod_tau.py -- production (sampling, HEURISTIC coverage; search/NOTES.md section 92) for ONE tau-alive real plug-in shape: the exact row-pattern classes (tau_patterns) give the allocation tokens n.m:n.m:n.m:n.m;
gs_samp.run is called with GS_PATTERNS=1 (tau-exact row_patterns3) and GS_EXACT_CAP (pattern-restricted exact lists when small: EXHAUSTED verdicts per (allocation, pattern)).
usage: prod_tau.py MODE Kmods r hg SECS SEED OUTPREFIX [LA LC [MAXTOK]]"""
import sys, os, json, time, itertools, collections
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault('GS_PATTERNS', '1')
os.environ.setdefault('GS_EXACT_CAP', '3e6')
from gs_samp import *
import tau_patterns as tp

def tokens_for(sh, maxtok=None):
    out, S = tp.patterns(sh)
    if out is None: return None
    res, Fu, inv, S = out
    byf = collections.defaultdict(list)
    for idx, f in enumerate(inv): byf[int(f)].append(idx)
    allocs = collections.Counter()
    for (a, b, c, d) in res:
        for fs in itertools.product(byf[a], byf[b], byf[c], byf[d]):
            al = tuple(sorted((int(S[i][0]), int(S[i][1])) for i in fs))
            allocs[al] += 1
    toks = sorted(allocs.items(), key=lambda x: (sum(n for n, m in x[0]), x[0]))
    return [al for al, cnt in toks]

if __name__ == '__main__':
    mode = sys.argv[1]; Kmods = json.loads(sys.argv[2]); r = int(sys.argv[3]); hg = json.loads(sys.argv[4])
    secs = float(sys.argv[5]); seed = int(sys.argv[6]); outp = sys.argv[7]
    LA = int(float(sys.argv[8])) if len(sys.argv) > 8 else 2000; LC = int(float(sys.argv[9])) if len(sys.argv) > 9 else 200000
    maxtok = int(sys.argv[10]) if len(sys.argv) > 10 else None
    sh = make_shape(mode, Kmods, r, hg, None)
    toks = tokens_for(sh)
    if toks is None: print('pattern enumeration too large'); sys.exit(1)
    print('tau-feasible allocation tokens:', len(toks), flush=True)
    if maxtok: toks = toks[:maxtok]
    txt = ';'.join(':'.join(f'{n}.{m}' for n, m in al) for al in toks)
    allocs = parse_allocs(txt, LA, LC)
    print(run(sh, allocs, secs, seed, outp, log=lambda *a: print(*a, flush=True), use_patterns=True, maxinst=int(float(os.environ.get('GS_MAXINST', '20000')))), flush=True)

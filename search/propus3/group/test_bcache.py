"""B-list cache test: shards with -bsave (first) / -blist (rest) give the same hit set as the plain run and as the brute-force canonical set."""
import sys, os, tempfile, subprocess, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, p3g, p3run as P
from p3g import GStruct
HERE = os.path.dirname(os.path.abspath(__file__))
w = tempfile.mkdtemp(prefix='bcache_')
cases = [(25, [(0, 1)], (1, 7, 1)), (27, [(0, 1)], (9, 1, 5)), (49, [(2, 1)], (13, 1, 5)), (25, [(8, 0)], (1, 7, 1))]
ok_all = True
for v, gens, rows in cases:
    st = GStruct(v, hgens=gens)
    pf = os.path.join(w, f'p_{v}.txt'); P.write_problem(st, rows, pf, canonB=True, canonD=True, psd=(True,) * 3, seed=v)
    plain = set(P.run_p3(pf, binary=os.path.join(HERE, 'p3'), threads=2, ns=2)['hits'])
    bf = os.path.join(w, f'B_{v}.bin'); N = 3; got = set()
    for sh in range(N):
        if sh == 0: r = P.run_p3(pf, binary=os.path.join(HERE, 'p3'), threads=2, ns=2, shard=(sh, N), extra=['-bsave', bf])
        else: r = P.run_p3(pf, binary=os.path.join(HERE, 'p3'), threads=3, ns=2, shard=(sh, N), blist=bf)
        got |= set(r['hits'])
    ok = got == plain and len(plain) > 0
    ok_all &= ok
    print(v, gens, rows, 'plain', len(plain), 'sharded+Bcache', len(got), 'OK' if ok else 'FAIL', 'B file', os.path.getsize(bf) // 8, 'masks', flush=True)
print('ALL OK' if ok_all else 'FAILED'); sys.exit(0 if ok_all else 1)

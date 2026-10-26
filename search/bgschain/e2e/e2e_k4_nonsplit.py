"""END-TO-END control, k = 4, NON-SPLIT G = Z4^2 (S = 2G, Q = Z2^2, n = 19, order 76): engine solution -> border -> certificate -> repository verify
(verify.build_checked).   usage: python3 e2e/e2e_k4_nonsplit.py [CERT.json]   (default results/k4_nonsplit_Z4x4_certificate.json under the current directory;
the engine solution record is written next to it as ..._solution.json)"""
import sys, os, json, tempfile, shutil
HERE = os.path.dirname(os.path.abspath(__file__)); ENG = os.path.join(HERE, '..')
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join('results', 'k4_nonsplit_Z4x4_certificate.json')
sys.path.insert(0, ENG); sys.path.insert(0, os.path.join(ENG, 'tests')); sys.path.insert(0, HERE)
import numpy as np, hlib, hchain, nonsplit_lib as NL, cert_lib as CL, border4
m2 = [4, 4]; S2 = [[2, 0], [0, 2]]; g2 = [NL.T2]
pats = NL.all_patterns(m2, S2, g2, -12, 4)
print('patterns', len(pats))
done = 0
for pt in pats:
    cfg = NL.real_cfg(m2, S2, g2, pt, -12, 4); cfg.update(levels=[{'ker': S2}], bdepth=1, sym_swap=True)
    d = tempfile.mkdtemp(); ch = hchain.Chain(cfg, d, verbose=False); data = ch.prepare()
    found = []
    for b in range(len(data['w'])):
        st, ctx = ch.run_branch(data, b); found += ctx['found']
    shutil.rmtree(d)
    if not found: continue
    for ms in found[:3]:
        rec = ch.sequences_json(ms)
        res = border4.find_border([np.array(y) for y in rec['ys']], k=4)
        if res is None: print('pattern', [(y[0], y[1]) for y in pt], 'no border found for this solution'); continue
        K, E, F = res
        cert = CL.make_certificate(rec, K, E, F, -12, 4)
        n, digest, cert = CL.verify_with_repo(cert, OUT)
        print(f'REPO VERIFY OK: order {n}, H H^T = {n} I, sha256 {digest}; group {cert["group"]}, subgroup {cert["subgroup_generators"]}, rowsums {cert["rowsums"]}')
        json.dump(rec, open(os.path.splitext(OUT)[0] + '_solution.json', 'w'))
        done = 1; break
    if done: break
print('E2E k4 non-split:', 'DONE' if done else 'NO SOLUTION/BORDER')

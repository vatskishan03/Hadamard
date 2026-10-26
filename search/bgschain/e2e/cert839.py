#!/usr/bin/env python3
"""For an n = 839 solution record (a line of solutions.jsonl): independent exact verification over G, border from z7search, certificate, and the
repository's own verification (verify.build_checked from ../../../verify/verify.py, in place) of H H^T = 3356 I.

    python3 e2e/cert839.py SOLUTIONS.jsonl [--line N] [--out CERT.json] [--solution-json PATH]

--out            certificate file written (repository format, with sha256; default results/certificate_3356.json under the current directory)
--solution-json  optional: a record of the solution (group, generators, sequences, row sums, pattern, sha256)
"""
import sys, os, json, argparse
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import cert_lib as CL, border8
ap = argparse.ArgumentParser()
ap.add_argument('solutions'); ap.add_argument('--line', type=int, default=0)
ap.add_argument('--out', default=os.path.join('results', 'certificate_3356.json')); ap.add_argument('--solution-json')
args = ap.parse_args()
mods = [4, 4, 4, 13]; v = 832; k = 8
rec = json.loads(open(args.solutions).read().strip().split('\n')[args.line])
seqs = [np.array([1 if c == '+' else -1 for c in s]) for s in rec['sequences']]
coords = np.array(np.unravel_index(np.arange(v), mods)).T
idx = lambda c: np.ravel_multi_index(tuple(np.moveaxis(c % mods, -1, 0)), mods)
# --- independent verification over G
h = np.array(rec['gens'][0]); img = idx(coords @ h.T)
assert all(np.array_equal(x[img], x) for x in seqs), 'not H-invariant'
Sset = (coords[:, 0] % 2 == 0) & (coords[:, 1] % 2 == 0) & (coords[:, 2] % 2 == 0)
cosets = (coords[:, 0] % 2) * 4 + (coords[:, 1] % 2) * 2 + (coords[:, 2] % 2)
ys = [[int(x[cosets == c].sum()) for c in range(8)] for x in seqs]
tot = np.zeros(v, dtype=np.int64)
for s in range(v):
    sh = idx(coords + coords[s]); tot[s] = sum(int(x @ x[sh]) for x in seqs)
target = np.where(Sset, -28, 4); target[0] = 4 * v
assert np.array_equal(tot, target), 'PAF condition fails'
A = [y[0] for y in ys]; a = [y[1] for y in ys]
assert all(y[1:] == [y[1]] * 7 for y in ys)
d = [A[i] - a[i] for i in range(4)]
print('verified over G: H-invariance, coset sums', ys[0], '..., sum PAF = -28 on S\\{0}, +4 off S;  pattern a', a, 'd', d, 'rowsums', [int(x.sum()) for x in seqs])
K, E, F = border8.make_border(a, d)
cert = CL.make_certificate(dict(order=3356, group=mods, subgroup_generators=rec['subgroup_generators'], gens=rec['gens'], sequences=rec['sequences']), K, E, F, -28, 4)
n, digest, cert = CL.verify_with_repo(cert, args.out)
print(f'REPO VERIFY OK (verify.build_checked): order {n}, H H^T = {n} I, sha256 {digest}')
print('wrote', args.out)
if args.solution_json:
    sol = dict(order=3356, n=839, group=mods, subgroup_generators=rec['subgroup_generators'], H_generator=rec['gens'][0], sequences=rec['sequences'], rowsums=[int(x.sum()) for x in seqs],
               pattern=dict(a=a, d=d, A=A, coset_sums=ys), autocorrelation=dict(subgroup=-28, elsewhere=4), sha256=digest, verified='over G (exact) and by verify.build_checked')
    os.makedirs(os.path.dirname(os.path.abspath(args.solution_json)), exist_ok=True)
    json.dump(sol, open(args.solution_json, 'w')); print('wrote', args.solution_json)

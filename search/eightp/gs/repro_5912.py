"""repro_5912.py -- deterministic replay of the hit of H(5912): round 33 of the production with base seed 1002 (worker w2), token (0, 30, 30, 60) (README, "Order 5912").
Regenerates the level-1 instances of that round (prod_5912.round_instances), runs level 2, compares the solution with the recorded K-level blocks, checks it at G level (exact integer PAF),
writes the certificate (construction spence, group [211, 7]) with gs_cert.certify_solution (H built by the verify.py of repo_copy, a copy of verify/verify.py; H H^T = 5912 I, sha256, then the CLI verify/verify.py) and compares its
sequences with REF.   Needs the builds and repo_copy of the README (Build).
usage: repro_5912.py [--fast] [--no-build] [REF]
  --fast      level 2 on the recorded instance only (level 1 of the round is still regenerated)
  --no-build  stop after the G-level check: write out/spence_5912_repro5912.json without sha256 (no matrix is built)
  REF         certificate to compare with (default ../../../hadamard/5912/certificate.json, if it exists)"""
import sys, os, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from gs_samp import make_shape, parse_allocs, process_instances
from gsk import verify_g
from gl1 import L1
from gs_cert import REPO, build_cert, certify_solution
from prod_5912 import SHAPE, TOKENS, LA, LC, round_instances

HIT = os.path.join(HERE, '..', 'data', 'gs', 'p5912', 'prod_w2.HIT_1793224543.json')
BASE, ROUND = 1002, 33


def main():
    fast = '--fast' in sys.argv; build = '--no-build' not in sys.argv
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    ref = args[0] if args else os.path.join(HERE, '..', '..', '..', 'hadamard', '5912', 'certificate.json')
    rec = json.load(open(HIT)); hit = rec['blocks'][0]                       # four (s, t, e) lists of length 211
    assert (rec['base'], rec['round']) == (BASE, ROUND)
    assert not build or os.path.isdir(os.path.join(REPO, 'verify')), 'repo_copy is missing: see the README (Build)'
    sh = make_shape(*SHAPE); l1 = L1(sh); allocs = parse_allocs(TOKENS, LA, LC)
    t0 = time.time()
    rng, al, cnt, S, E = round_instances(sh, l1, allocs, {}, BASE, ROUND)
    print(f'round {ROUND} of base seed {BASE}, token {al}: {cnt} level-1 instances, {len(S)} distinct ({time.time()-t0:.0f} s)', flush=True)
    if fast:
        i = next((i for i in range(len(S)) if all(np.array_equal(S[i, b], hit[b][0]) and np.array_equal(E[i, b], hit[b][2]) for b in range(4))), None)
        assert i is not None, 'recorded instance not regenerated: level 1 of round 33 differs from the production'
        S, E = S[i:i + 1], E[i:i + 1]; print(f'recorded instance: number {i}', flush=True)
    ncomp, sols, stat, tl = process_instances(sh, S, E, rng, tl2=30.0)
    print(f'level 2: {len(S)} instances, {ncomp} completions (one solution, doubled by t -> -t), {len(sols)} solutions kept ({tl:.0f} s)', flush=True)
    same = [all(np.array_equal(bl[b][i], hit[b][i]) for b in range(4) for i in range(3)) for bl in sols]
    assert any(same), 'the replay does not reproduce the recorded solution'
    bl = sols[same.index(True)]
    ok, gseqs = verify_g(sh, bl); assert ok, 'G-level check failed'
    print(f'solution: identical (s, t, e) blocks to the recorded hit; G-level exact check: sum PAF = 5908 at 0, -4 elsewhere, row sums {tuple(int(np.sum(x)) for x in gseqs)}', flush=True)
    if not build:
        cert = build_cert(sh, bl); os.makedirs(os.path.join(HERE, 'out'), exist_ok=True)
        json.dump(cert, open(os.path.join(HERE, 'out', 'spence_5912_repro5912.json'), 'w')); print('certificate without sha256: out/spence_5912_repro5912.json'); res = None
    else:
        res = certify_solution(sh, bl, tag='repro5912')
        print('certificate', res, flush=True); assert res['verify_ok'], 'verify.py failed'
        cert = json.load(open(os.path.join(res['dir'], 'certificate.json')))
    if os.path.exists(ref):
        r = json.load(open(ref)); keys = ('order', 'construction', 'group', 'sequences') + (('sha256',) if build else ())
        print(f'compared with {os.path.relpath(ref, HERE)} ({", ".join(keys)}): identical =', all(r[k] == cert[k] for k in keys), flush=True)
    else: print('no reference certificate at', os.path.relpath(ref, HERE))
    print(f'total {time.time()-t0:.0f} s')


if __name__ == '__main__': main()

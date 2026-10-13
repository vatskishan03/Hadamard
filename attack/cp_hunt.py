"""Estimate-then-search loop over cyclotomic-product structures (cp_sweep.structures) with the low-energy k-tree of
cp_lowe.py.

For each target v (skew: odd v of the skew census; plain: n of FRONTIER.md for Goethals-Seidel order 4n) and each
structure, three stage-1 classes are sampled (cp_lowe.stage1_sample, seed 1) and the Gaussian estimate of the number of
solutions of the first one inside low-energy lists of L members per non-skew block is computed. When it is at least
MINE and the predicted k-tree work (pairs AB + CD per solution) is at most 10^MAXWORK, cp_lowe.py is run on that class
(as a subprocess with a time limit). A v is left as soon as one structure gives a verified solution; skew hits are then
written to skew/<4v>/ by skew_certify.py.

  python3 attack/cp_hunt.py skew|plain v1,v2,.. [--L 1e7] [--mine 0.5] [--maxwork 11.5] [--timeout 3600] [--threads 3]
      -> $HADAMARD_WORK/cp_hunt.jsonl (one line per structure), logs in $HADAMARD_WORK/logs/cp_hunt_<v>_*.log
"""
import argparse, json, math, os, subprocess, sys, time
import numpy as np
from pathlib import Path as _P
HERE = _P(__file__).resolve().parent
sys.path[:0] = [str(HERE)]
import workdir as W
import cp_lowe as CL
import cp_sweep as SW


def solved(v, skew):
    path = W.hits() if skew else W.hits('had_hits.jsonl')
    if not os.path.exists(path): return False
    for l in open(path):
        try: r = json.loads(l)
        except Exception: continue
        if r.get('n') == v and r.get('skew', False) == skew and r.get('hadamard'): return True
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=('skew', 'plain')); ap.add_argument('vs')
    ap.add_argument('--L', type=float, default=1e7); ap.add_argument('--mine', type=float, default=0.5)
    ap.add_argument('--maxwork', type=float, default=11.5); ap.add_argument('--timeout', type=float, default=3600)
    ap.add_argument('--threads', type=int, default=3); ap.add_argument('--samples', type=int, default=60000)
    ap.add_argument('--minm', type=float, default=1e5); ap.add_argument('--maxm', type=float, default=3e10)
    ap.add_argument('--nsig', type=int, default=3000)
    a = ap.parse_args()
    skew = a.mode == 'skew'
    out = W.WORK / 'cp_hunt.jsonl'
    logd = W.WORK / 'logs'; logd.mkdir(exist_ok=True)
    done = set()
    if out.exists():
        for l in open(out):
            r = json.loads(l); done.add((r['v'], r['mode'], r['q'], r['e'], r['k'], r['twist']))
    for v in [int(x) for x in a.vs.split(',')]:
        if solved(v, skew): continue
        for q, m, e, F, k, t, S, mem in SW.structures(v, skew, a.minm, a.maxm):
            if solved(v, skew): break
            if (v, a.mode, q, e, k, t) in done: continue
            base = dict(v=v, mode=a.mode, q=q, m=m, e=e, F=F, k=k, K=len(S.K), twist=t, reps=len(S.kreps),
                        mem0=round(math.log10(mem), 1))
            t0 = time.time()
            if len(S.kreps) > 22:                                  # 2^reps sign patterns per spike placement in stage 1
                continue
            try:
                O = CL.Orbits(S, t)
                quads = CL.stage1_sample(S, 3, np.random.default_rng(1), O=O)
                if not quads:
                    rec = dict(base, classes=0)
                else:
                    P = CL.prepare(S, O, quads[0], a.samples, np.random.default_rng(1))
                    if S.skew and len(P['Pm0']) == 0: raise RuntimeError('no skew block member')
                    est, rank, off = CL.estimate_in_lists(P, lambda kk: a.L)
                    sz = sorted(P['sizes']); Lb = [min(a.L, x) for x in sz[1:]]
                    work = math.log10(sz[0] * Lb[0] + Lb[1] * Lb[2]) - est
                    rec = dict(base, classes=len(quads), y=quads[0], rows=[int(sum(yy)) for yy in quads[0]], logE=round(P['est0'], 2),
                               logE_lists=round(est, 2), rank=rank, sizes=[round(math.log10(max(1, x)), 1) for x in P['sizes']],
                               logwork=round(work, 2))
                    if est >= a.mine and work <= a.maxwork:
                        log = logd / f'cp_hunt_{v}_{q}_{e}_{k}_{t}.log'
                        cmd = [sys.executable, '-u', str(HERE / 'cp_lowe.py'), str(v), '--q', str(q), '--e', str(e), '--kgen', str(k),
                               '--twist', str(t), '--sample-stage1', '3', '--config', '0', '--L', f'{a.L:g}', '--nsig', str(a.nsig),
                               '--threads', str(a.threads), '--K', '3'] + ([] if skew else ['--plain'])
                        try:
                            with open(log, 'w') as fh:
                                subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, timeout=a.timeout)
                            status = 'done'
                        except subprocess.TimeoutExpired:
                            status = 'timeout'
                        txt = open(log).read()
                        nver = 0
                        for line in txt.splitlines():
                            if line.endswith('verified solutions') or ' verified solutions [' in line:
                                nver = int(line.split()[0])
                        rec.update(search=status, cmd=' '.join(cmd[1:]), verified=nver, log=str(log))
                        if nver == 0 and 'hit: autocorrelations ok' in txt and solved(v, skew): rec['verified'] = 1
                        if rec['verified'] and skew:
                            subprocess.run([sys.executable, str(HERE / 'skew_certify.py'), str(v), '0'], timeout=1800)
            except (RuntimeError, MemoryError, AssertionError, ValueError, StopIteration) as ex:
                rec = dict(base, error=repr(ex)[:100])
            rec['sec'] = round(time.time() - t0)
            print(json.dumps(rec), flush=True)
            with open(out, 'a') as f: f.write(json.dumps(rec) + '\n')


if __name__ == '__main__':
    main()

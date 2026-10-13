"""Breadth over stage-1 classes: run cp_lowe.py on classes c0, c0+1, ... of one structure with a small number of split
values each (the best-ranked split values of a class hold the densest share of its solutions), and stop at the first
verified solution.

  python3 attack/cp_breadth.py v --q Q --e E [--kgen k] [--twist t] [--field p,k] [--sample-stage1 N]
                               [--classes 0-35] [--L 2e7,1e7,1e7] [--nsig 300] [--threads 3] [--timeout 7200]
      logs in $HADAMARD_WORK/logs/cp_breadth_<v>_<class>.log; skew hits are certified (skew_certify[_group].py).
"""
import argparse, os, subprocess, sys, time
from pathlib import Path as _P
HERE = _P(__file__).resolve().parent
sys.path[:0] = [str(HERE)]
import workdir as W


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('v', type=int); ap.add_argument('--q', type=int, required=True); ap.add_argument('--e', type=int, required=True)
    ap.add_argument('--kgen', type=int, default=1); ap.add_argument('--twist', type=int, default=0)
    ap.add_argument('--field', default=''); ap.add_argument('--sample-stage1', type=int, default=0)
    ap.add_argument('--classes', default='0-9'); ap.add_argument('--L', default='2e7,1e7,1e7')
    ap.add_argument('--nsig', type=int, default=300); ap.add_argument('--threads', type=int, default=3)
    ap.add_argument('--timeout', type=float, default=7200); ap.add_argument('--plain', action='store_true')
    a = ap.parse_args()
    lo, hi = (int(x) for x in a.classes.split('-'))
    logd = W.WORK / 'logs'; logd.mkdir(exist_ok=True)
    for c in range(lo, hi + 1):
        log = logd / f'cp_breadth_{a.v}_{c}.log'
        cmd = [sys.executable, '-u', str(HERE / 'cp_lowe.py'), str(a.v), '--q', str(a.q), '--e', str(a.e), '--kgen', str(a.kgen),
               '--twist', str(a.twist), '--config', str(c), '--L', a.L, '--nsig', str(a.nsig), '--threads', str(a.threads), '--K', '1']
        if a.field: cmd += ['--field', a.field]
        if a.sample_stage1: cmd += ['--sample-stage1', str(a.sample_stage1)]
        if a.plain: cmd += ['--plain']
        t0 = time.time()
        try:
            with open(log, 'w') as fh: subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, timeout=a.timeout)
            status = 'done'
        except subprocess.TimeoutExpired:
            status = 'timeout'
        txt = open(log).read()
        nver = sum(1 for l in txt.splitlines() if l.startswith('hit: invariant') and ': True [' in l)
        est = [l for l in txt.splitlines() if 'inside the kept lists' in l]
        print(f'class {c}: {status}, {nver} verified, {est[0].split(":", 1)[1].strip() if est else "no estimate"} [{time.time() - t0:.0f}s]', flush=True)
        if 'IndexError' in txt or 'list index out of range' in txt: break                  # past the last class
        if nver:
            if not a.plain:
                certifier = 'skew_certify_group.py' if a.field else 'skew_certify.py'
                subprocess.run([sys.executable, str(HERE / certifier), str(a.v), '0'], timeout=1800)
            break


if __name__ == '__main__':
    main()

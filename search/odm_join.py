"""Large-list version of odm_engine.py: the c parts of each H-invariant template are joined by join4x.c instead of
odmx's in-memory pair hash.

For every template printed by `odm_engine.py --dry` (b parts of weight m, row-sum screen passed): `odmx ... DUMP`
writes the four lists of c parts of each row-sum vector (vectors of autocorrelations and bc terms at the pair-orbit
representatives) in join4x's input format; join4x finds every quadruple summing to zero (exactly, with
memory-bounded partitions); each solution is decoded back to the four blocks with `odmx ... --decode`.

Usage: python3 search/odm_join.py t gens m [--bsub f] [--work DIR] [--k K] [--maxtpl N] [--first]
"""
import argparse
import os
import struct
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('t', type=int)
    ap.add_argument('gens')
    ap.add_argument('m', type=int)
    ap.add_argument('--bsub', type=int, default=0)
    ap.add_argument('--work', default=os.path.join(os.environ.get('HADAMARD_WORK', '/tmp/hadamard-work'), 'tmp'))
    ap.add_argument('--odmx', default='/tmp/hadamard-work/bin/odmx')
    ap.add_argument('--join4x', default='/tmp/hadamard-work/bin/join4x')
    ap.add_argument('--k', type=int, default=2, help='solutions per join')
    ap.add_argument('--maxtpl', type=int, default=10 ** 6)
    ap.add_argument('--first', action='store_true', help='stop at the first solution')
    ap.add_argument('--bjoin', action='store_true', help='templates via join4x (odm_engine --bjoin)')
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)
    cmd = [sys.executable, os.path.join(HERE, 'odm_engine.py'), str(a.t), a.gens, str(a.m), '--dry']
    if a.bsub:
        cmd += ['--bsub', str(a.bsub)]
    if a.bjoin:
        cmd += ['--bjoin']
    r = subprocess.run(cmd, capture_output=True, text=True)
    print(r.stderr.strip(), file=sys.stderr)
    tpls = [ln.split(' ', 1)[0] for ln in r.stdout.splitlines() if ln.strip()]
    print(f'{len(tpls)} templates', file=sys.stderr)
    total = 0
    for k, tpl in enumerate(tpls[:a.maxtpl]):
        prefix = os.path.join(a.work, f'odmj_{a.t}_{a.m}_{os.getpid()}_{k}')
        d = subprocess.run([a.odmx, str(a.t), a.gens, tpl, '1', prefix], capture_output=True, text=True)
        dumps = [ln.split()[1] for ln in d.stdout.splitlines() if ln.startswith('DUMP')]
        nsol = 0
        for base in dumps:
            j = subprocess.run([a.join4x, base + '.bin', str(a.k)], capture_output=True, text=True)
            sols = [list(map(int, ln.split()[1:5])) for ln in j.stdout.splitlines() if ln.startswith('SOL')]
            tot = [ln for ln in j.stdout.splitlines() if ln.startswith('TOTAL')]
            if sols:
                with open(base + '.idx', 'rb') as f:
                    n = struct.unpack('4i', f.read(16))
                    masks = [struct.unpack(f'{c}i', f.read(4 * c)) for c in n]
                for s in sols:
                    ms = [str(masks[i][s[i]]) for i in range(4)]
                    out = subprocess.run([a.odmx, str(a.t), a.gens, tpl, '1', '--decode'] + ms,
                                         capture_output=True, text=True).stdout
                    print(f'# template {tpl}')
                    print(out, end='')
                    sys.stdout.flush()
                    nsol += 1
            print(f'template {k}: {os.path.basename(base)} {tot[0] if tot else "?"}', file=sys.stderr)
            for ext in ('.bin', '.idx'):
                try:
                    os.remove(base + ext)
                except OSError:
                    pass
        total += nsol
        if nsol and a.first:
            break
    print(f'total solutions {total}', file=sys.stderr)


if __name__ == '__main__':
    main()

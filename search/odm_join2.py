"""odm_join.py with cached c lists: every list of H-invariant c parts is built once per (block type, b part, row
sum) by `odmc.c` (depth-first search with incremental vectors) and reused by all templates; each template's
row-sum vectors are then joined by join4x.

Usage: python3 search/odm_join2.py t gens m [--bsub f] [--bjoin] [--k K] [--first] [--maxtpl N]
"""
import argparse
import ast
import os
import struct
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import workdir as W  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('t', type=int)
    ap.add_argument('gens')
    ap.add_argument('m', type=int)
    ap.add_argument('--bsub', type=int, default=0)
    ap.add_argument('--bjoin', action='store_true')
    ap.add_argument('--k', type=int, default=2)
    ap.add_argument('--first', action='store_true')
    ap.add_argument('--maxtpl', type=int, default=10 ** 7)
    ap.add_argument('--work', default=W.tmpdir())
    ap.add_argument('--odmc', default=None,
                    help='binary built from odmc.c (default: compiled into $HADAMARD_WORK/bin on first use)')
    ap.add_argument('--join4x', default=None,
                    help='binary built from join4x.c (default: compiled into $HADAMARD_WORK/bin on first use)')
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)
    cmd = [sys.executable, os.path.join(HERE, 'odm_engine.py'), str(a.t), a.gens, str(a.m), '--dry']
    if a.bsub:
        cmd += ['--bsub', str(a.bsub)]
    if a.bjoin:
        cmd += ['--bjoin']
        if a.join4x:
            cmd += ['--join4x', a.join4x]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        sys.exit(r.stderr.strip() or f'odm_engine.py failed ({r.returncode})')
    print(r.stderr.strip(), file=sys.stderr)
    tpls = []
    for ln in r.stdout.splitlines():
        if not ln.strip():
            continue
        tstr, feas = ln.split(' ', 1)
        tpls.append((tstr, ast.literal_eval(feas)))
    print(f'{len(tpls)} templates', file=sys.stderr)
    cache = {}
    tag = f'{a.t}_{a.m}_{a.gens.replace(",", "_")}_{os.getpid()}'

    def block_entries(tstr, i):
        return [e for e in tstr.split(',') if e.split(':')[0] == str(i)]

    def get_list(tstr, i, rs):
        ents = block_entries(tstr, i)
        key = (1 if i == 1 else 2, tuple(sorted(e.split(':', 1)[1] for e in ents)), rs)
        if key in cache:
            return cache[key]
        fn = os.path.join(a.work, f'odmc_{tag}_{len(cache)}.lst')
        relab = ','.join(f'{i}:' + e.split(':', 1)[1] for e in ents) or f'{i}:0:0'
        subprocess.run([(a.odmc or W.bin('odmc')), str(a.t), a.gens, str(i), relab, str(rs), fn], check=True, capture_output=True)
        with open(fn, 'rb') as f:
            n, D = struct.unpack('2i', f.read(8))
            vec = f.read(2 * n * D)
            masks = struct.unpack(f'{n}Q', f.read(8 * n)) if n else ()
        os.remove(fn)
        cache[key] = (n, D, vec, masks)
        return cache[key]

    total = 0
    for k, (tstr, feas) in enumerate(tpls[:a.maxtpl]):
        nsol = 0
        for (r2, r3, r4) in feas:
            L = [get_list(tstr, 1, 0), get_list(tstr, 2, r2), get_list(tstr, 3, r3), get_list(tstr, 4, r4)]
            if any(x[0] == 0 for x in L):
                continue
            D = L[0][1]
            fn = os.path.join(a.work, f'odmj2_{tag}.bin')
            with open(fn, 'wb') as f:
                f.write(struct.pack('5i', L[0][0], L[1][0], L[2][0], L[3][0], D))
                f.write(b'\0' * (2 * D))
                for x in L:
                    f.write(x[2])
            j = subprocess.run([(a.join4x or W.bin('join4x')), fn, str(a.k)], capture_output=True, text=True)
            os.remove(fn)
            sols = [list(map(int, ln.split()[1:5])) for ln in j.stdout.splitlines() if ln.startswith('SOL')]
            tot = [ln for ln in j.stdout.splitlines() if ln.startswith('TOTAL')]
            print(f'template {k} rows (0, {r2}, {r3}, {r4}): lists {[x[0] for x in L]} {tot[0] if tot else "?"}',
                  file=sys.stderr)
            for s in sols:
                rows = []
                for i in range(4):
                    ents = block_entries(tstr, i + 1)
                    relab = ','.join(ents) or f'{i + 1}:0:0'
                    out = subprocess.run([(a.odmc or W.bin('odmc')), str(a.t), a.gens, str(i + 1), relab, '0', '--decode',
                                          str(L[i][3][s[i]])], capture_output=True, text=True).stdout.strip()
                    rows.append(out)
                print(f'# template {tstr}')
                print(f'SOL t={a.t} m={a.m} rows (0, {r2}, {r3}, {r4})')
                print('\n'.join(rows))
                sys.stdout.flush()
                nsol += 1
        total += nsol
        if nsol and a.first:
            break
    print(f'total solutions {total}; {len(cache)} distinct c lists', file=sys.stderr)


if __name__ == '__main__':
    main()

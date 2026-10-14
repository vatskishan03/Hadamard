"""Structure files for orbit_ls.c (tabu search on multiplier-orbit structures) and a driver that checks the hits.

  python3 search/orbit_ls.py v gens family [--seed S] [--maxsec T] [--nsol N] [--R R]
family: skew (k1,f1,f1,f1; skew Goethals-Seidel), propus (s1,f2,f1), gs (f1,f1,f1,f1), williamson (s1,s1,s1,s1),
lp (f1,f1, T = -2), dopt (f1,f1, T = 2). gens: comma-separated generators of H (1 for no symmetry).
Hits are checked (sum of weighted PAFs) and appended to $HADAMARD_WORK/orbit_ls_hits.jsonl.
"""
import argparse, json, os, subprocess, sys, tempfile
import numpy as np
from pathlib import Path as _P
sys.path[:0] = [str(_P(__file__).resolve().parent)]
import workdir as W

FAMILIES = {'skew': ('k1,f1,f1,f1', 0), 'propus': ('s1,f2,f1', 0), 'gs': ('f1,f1,f1,f1', 0),
            'williamson': ('s1,s1,s1,s1', 0), 'lp': ('f1,f1', -2), 'dopt': ('f1,f1', 2)}


def subgroup(v, gens):
    H = {1}
    while True:
        H2 = H | {h * g % v for h in H for g in gens}
        if H2 == H: return sorted(H)
        H = H2


def structure_file(v, gens, family, path):
    types, T = FAMILIES[family]
    blocks = types.split(',')
    H = subgroup(v, gens)
    oid = [-1] * v; orbs = []
    for x in range(v):
        if oid[x] < 0:
            o = sorted({x * h % v for h in H})
            for y in o: oid[y] = len(orbs)
            orbs.append(o)
    mir = [oid[(-o[0]) % v] for o in orbs]
    Hm = sorted(set(H) | {(-h) % v for h in H})
    seen = [False] * v; classes = []
    for s in range(1, v):
        if not seen[s]:
            cl = {s * h % v for h in Hm}
            for y in cl: seen[y] = True
            classes.append((len(cl), s))
    with open(path, 'w') as f:
        f.write(f"{v} {len(blocks)} {T}\n")
        for b in blocks: f.write(f"{b[0]} {b[1]}\n")
        f.write(f"{len(orbs)}\n")
        for o in orbs: f.write(f"{len(o)} " + ' '.join(map(str, o)) + "\n")
        f.write("mirror " + ' '.join(map(str, mir)) + "\n")
        f.write(f"{len(classes)}\n")
        for n, s in classes: f.write(f"{n} {s}\n")
    return H, orbs, classes, blocks, T


def paf(x):
    v = len(x)
    return np.array([int(np.dot(x, np.roll(x, -s))) for s in range(v)])


def check(seqs, blocks, T):
    tot = sum(int(b[1]) * paf(x) for x, b in zip(seqs, blocks))
    return bool(np.all(tot[1:] == T)), int(tot[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('v', type=int); ap.add_argument('gens'); ap.add_argument('family', choices=sorted(FAMILIES))
    ap.add_argument('--seed', type=int, default=1); ap.add_argument('--maxsec', type=float, default=60)
    ap.add_argument('--nsol', type=int, default=1); ap.add_argument('--R', type=int, default=20000)
    ap.add_argument('--TT', type=int, default=0)
    a = ap.parse_args()
    gens = [int(g) for g in a.gens.split(',')]
    fd, path = tempfile.mkstemp(suffix='.txt', dir=W.tmpdir()); os.close(fd)
    try:
        H, orbs, classes, blocks, T = structure_file(a.v, gens, a.family, path)
        print(f"v={a.v} H={gens} |H|={len(H)} orbits {len(orbs)} classes {len(classes)} family {a.family}", flush=True)
        p = subprocess.Popen([W.bin('orbit_ls'), path, str(a.seed), str(a.maxsec), str(a.R), str(a.TT), str(a.nsol)],
                             stdout=subprocess.PIPE, text=True)
        for line in p.stdout:
            line = line.strip()
            if line.startswith('SOL'):
                parts = line.split()
                seqs = [np.array([1 if c == '+' else -1 for c in s]) for s in parts[1:1 + len(blocks)]]
                ok, at0 = check(seqs, blocks, T)
                print(f"hit: {a.family} v={a.v} H={gens} verified {ok} ({' '.join(parts[1 + len(blocks):])})", flush=True)
                if ok:
                    rec = dict(v=a.v, gens=gens, family=a.family, blocks=blocks, T=T, seqs=parts[1:1 + len(blocks)],
                               rows=[int(x.sum()) for x in seqs])
                    with open(W.hits('orbit_ls_hits.jsonl'), 'a') as f: f.write(json.dumps(rec) + '\n')
            else:
                print(line, flush=True)
        p.wait()
    finally:
        os.remove(path)


if __name__ == '__main__':
    main()

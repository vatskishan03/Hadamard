"""validate the symmetry reductions (uA, rB, sym, s3) of qp37r by orbit coverage: python3 symtest.py SHAPE [extra flags]
full run: all normalised solutions; reduced run: flags uA=1 rB=1 sym=1 s3=1; check that the G-orbits of the reduced solutions cover the full set exactly (G = units x reflections of B,C,D x S3(B,C,D))."""
import sys, os, subprocess, itertools
fn = sys.argv[1]
def run(flags):
    out = subprocess.run([os.environ.get('QP', './qp37r'), fn, 'real=1', 'l2c=1', 'maxprint=10000000'] + flags, capture_output=True, text=True).stdout
    sols = [tuple(l.split()[1:]) for l in out.split('\n') if l.startswith('SOL')]
    cnt = [l for l in out.split('\n') if l.startswith('COUNT')]
    return sols, cnt
lines = open(fn).read().split('\n')
nu, q, K = map(int, lines[0].split()); neg = list(map(int, lines[2].split()))
# unit perms at the tail
toks = open(fn).read().split(); pos = 3 + 2 * nu + nu ** 3 + 2 * nu * nu
nU = int(toks[pos]); U = [list(map(int, toks[pos + 1 + i * nu: pos + 1 + (i + 1) * nu])) for i in range(nU)]
U = [list(range(nu))] + U
def perm_block(b, p):   # new[o] = old[p[o]]
    return '|'.join(''.join(s[p[o]] for o in range(nu)) for s in b.split('|'))
def refl_block(b):
    x0, x1, x2 = b.split('|')
    return '|'.join([''.join(x0[neg[o]] for o in range(nu)), ''.join(x2[neg[o]] for o in range(nu)), ''.join(x1[neg[o]] for o in range(nu))])
def orbit(sol):
    res = set()
    for p in U:
        A, B, C, D = [perm_block(b, p) for b in sol]
        for rb, rc, rd in itertools.product((0, 1), repeat=3):
            B2 = refl_block(B) if rb else B; C2 = refl_block(C) if rc else C; D2 = refl_block(D) if rd else D
            for perm in itertools.permutations((B2, C2, D2)): res.add((A,) + perm)
    return res
full, c1 = run([]); red, c2 = run(['uA=1', 'rB=1', 'sym=1', 's3=1'])
print('full', len(full), c1, ' reduced', len(red), c2)
fs = set(full); cover = set()
for s in red:
    assert s in fs, 'reduced solution not in full set'
    cover |= orbit(s)
print('orbit closure of reduced =', len(cover), ' full =', len(fs), ' equal:', cover == fs, ' closure subset of full:', cover <= fs)

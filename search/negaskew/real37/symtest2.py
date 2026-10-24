"""validate tF=1 (translations by H_K-fixed points, reflection, sign) by orbit coverage: python3 symtest2.py SHAPE"""
import sys, os, subprocess, itertools
fn = sys.argv[1]
def run(flags):
    out = subprocess.run([os.environ.get('QP', './qp37r'), fn, 'real=1', 'l2c=1', 'maxprint=10000000'] + flags, capture_output=True, text=True).stdout
    sols = [tuple(l.split()[1:]) for l in out.split('\n') if l.startswith('SOL')]
    cnt = [l for l in out.split('\n') if l.startswith('COUNT')]
    return sols, cnt
toks = open(fn).read().split(); nu = int(toks[1 - 1]); q = int(toks[1])
nu = int(toks[0]); lines = open(fn).read().split('\n'); neg = list(map(int, lines[2].split()))
pos = 3 + 2 * nu + nu ** 3 + 2 * nu * nu
nU = int(toks[pos]); U = [[*range(nu)]] + [list(map(int, toks[pos + 1 + i * nu: pos + 1 + (i + 1) * nu])) for i in range(nU)]
pos2 = pos + 1 + nU * nu
nG = int(toks[pos2]); G = [(0, list(range(nu)))] + [(int(toks[pos2 + 1 + i * (nu + 1)]), list(map(int, toks[pos2 + 2 + i * (nu + 1): pos2 + 1 + (i + 1) * (nu + 1)]))) for i in range(nG)]
def flip(s): return ''.join('2' if c == '0' else '0' for c in s)
def norm(b):
    x0, x1, x2 = b.split('|')
    if x0[0] == '2': x0, x1, x2 = flip(x0), flip(x1), flip(x2)
    return '|'.join([x0, x1, x2])
def perm_block(b, p):
    return '|'.join(''.join(s[p[o]] for o in range(nu)) for s in b.split('|'))
def gblock(b, refl, p):
    x0, x1, x2 = b.split('|')
    if refl: x1, x2 = x2, x1
    return norm('|'.join(''.join(s[p[o]] for o in range(nu)) for s in (x0, x1, x2)))
def orbit(sol):
    res = set()
    for u in U:
        A = perm_block(sol[0], u); BCD = [perm_block(b, u) for b in sol[1:]]
        imgs = []
        for b in BCD:
            imgs.append({gblock(b, refl, p) for refl, p in G})
        for B2 in imgs[0]:
            for C2 in imgs[1]:
                for D2 in imgs[2]:
                    for perm in itertools.permutations((B2, C2, D2)): res.add((A,) + perm)
    return res
full, c1 = run(['uA=0'] if False else []); red, c2 = run(['uA=1', 'tF=1', 'sym=1', 's3=1'])
print('full', len(full), c1, ' reduced', len(red), c2)
fs = set(full)
cover = set()
for s in red:
    s2 = (s[0],) + tuple(norm(b) for b in s[1:])
    assert s2 in fs, 'reduced solution (normalised) not in full set'
    cover |= orbit(s2)
print('orbit closure of reduced =', len(cover), ' full =', len(fs), ' equal:', cover == fs, ' closure subset of full:', cover <= fs)

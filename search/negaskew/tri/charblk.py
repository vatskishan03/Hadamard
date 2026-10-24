"""Character-type blocks (search/NOTES.md section 92: none beyond chance): count the blocks of a solution set whose coset exponents are linear modulo 4 in the coset index, c_j = i^(a j + b) (a, b mod 4; a = 2 contains the Legendre-type blocks c_j = +-(-1)^j, a = 0 the constants), per position and per kind.  The coset index j is that of the tools of this directory (j -> g^j H, g the least primitive root).
usage: python3 charblk.py v k FILE        (FILE: SOL lines of ../qskew or of tri2)"""
import sys, os, collections; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trilib import *
v = int(sys.argv[1]); k = int(sys.argv[2]); d, sols = load_sols(sys.argv[3], v, k); e = d['e']
def lin(b):
    for a in range(4):
        for c in range(4):
            if all(b[1][j] == (a * j + c) % 4 for j in range(e)): return (a, c)
    return None
cnt = collections.Counter(); tot = collections.Counter()
for s in sols:
    for X in range(4):
        t = lin(s[X]); tot[X] += 1
        if t is not None and (X > 0 or True): cnt[(X, t)] += 1
print('v', v, 'k', k, 'e', e, 'solutions', len(sols), '; blocks with c_j = i^(aj+b), per position (0 = A, 1..3 = B, C, D):', dict(sorted(cnt.items())) if cnt else 'none')

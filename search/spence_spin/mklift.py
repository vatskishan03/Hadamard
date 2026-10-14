"""Write the lift problem file for liftsolve.c.  usage: mklift.py <wit.npy> <out.txt> [planted_seed]
planted_seed: replace the target by the keys of random signs s* (sA[0]=+1) and print s* (completeness test)."""
import sys, numpy as np
from lift_poly import Poly, natural_linear, build, make
S = build(3); LQ = make(S, 'Z2^2xZ13')
w = np.load(sys.argv[1]); a, b = w[:LQ.RK], w[LQ.RK:]
P = Poly(S, LQ, a, b)
T = np.full(S.D, -4, dtype=np.int64)
planted = None
if len(sys.argv) > 3:
    rng = np.random.default_rng(int(sys.argv[3]))
    sA = rng.choice([-1, 1], P.nA); sA[0] = 1; sB = rng.choice([-1, 1], P.nB)
    T = P.keys_poly(sA, sB); planted = (sA, sB)
    ke, _, _ = P.keys_exact(sA, sB); assert np.array_equal(ke, T)
E, rhs, ok, info = natural_linear(P, T)
assert ok, info
if planted is not None: assert np.array_equal(E @ np.concatenate(planted), rhs)
nA, nB = P.nA, P.nB
with open(sys.argv[2], 'w') as f:
    f.write(f"{nA} {nB} {len(E)} {S.D} 1\n")
    for e in range(len(E)): f.write(' '.join(map(str, E[e, :nA])) + '\n')
    for e in range(len(E)): f.write(' '.join(map(str, E[e, nA:])) + '\n')
    f.write(' '.join(map(str, rhs)) + '\n')
    for arr in (T, P.cA, P.cB): f.write(' '.join(map(str, arr)) + '\n')
    for c in range(S.D): f.write(' '.join(map(str, P.linA[c])) + '\n')
    for c in range(S.D): f.write(' '.join(map(str, P.linB[c])) + '\n')
    for c in range(S.D):
        for i in range(nA): f.write(' '.join(map(str, P.QA[c][i])) + '\n')
    for c in range(S.D):
        for i in range(nB): f.write(' '.join(map(str, P.QB[c][i])) + '\n')
print('nA', nA, 'nB', nB, 'E', len(E), info)
if planted is not None:
    mA = sum(1 << i for i in range(nA) if planted[0][i] == -1); mB = sum(1 << i for i in range(nB) if planted[1][i] == -1)
    print('PLANTED', mA, mB)

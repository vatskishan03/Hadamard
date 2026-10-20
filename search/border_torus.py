"""Bordered propus quadruples (NOTES.md §27) over the norm-one torus T of GF(q²), a cyclic group of order q + 1, for
v = q + 2 (q an odd prime; NOTES.md §63). A (symmetric), B, D are ±1 functions on T with
|Â|² + 2|B̂|² + |D̂|² = 4(q + 2) at every nontrivial character and row sums with a² + 2b² + d² = 4. The candidate
functions are χ(Tr t + c), χ(Tr(z t) + c) (χ the quadratic character of GF(q), z primitive in GF(q²)) and
η(t − β) (η the quadratic character of GF(q²), β ∈ λ·{1, z}), with both values at a zero. Every χ(Tr(αt) + c) is
± a shift of one of these, and η(t − β) = χ(1 + N(β) − Tr(β̄t)), so B and D need no other parameters (the spectrum is
shift-invariant); A must be symmetric, so every shift of every listed function is tried for A.

    python3 search/border_torus.py q     # distinct triples only at q = 3, 5; none for primes 7 ≤ q ≤ 47
"""
import sys, itertools
import numpy as np
from collections import defaultdict

q = int(sys.argv[1])
# GF(q^2) = GF(q)[s], s^2 = n (non-residue); elements (a, b) = a + b s
n = next(x for x in range(2, q) if pow(x, (q - 1) // 2, q) == q - 1)
def mul(x, y):
    return ((x[0] * y[0] + n * x[1] * y[1]) % q, (x[0] * y[1] + x[1] * y[0]) % q)
def pw(x, e):
    r = (1, 0)
    while e:
        if e & 1: r = mul(r, x)
        x = mul(x, x); e >>= 1
    return r
def chi(a):
    a %= q
    if a == 0: return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1
def norm(x): return (x[0] * x[0] - n * x[1] * x[1]) % q
def eta(x):  # quadratic character of GF(q^2) = chi(norm)
    return chi(norm(x))
# primitive element of GF(q^2)
order = q * q - 1
fac = [p for p in range(2, order + 1) if order % p == 0 and all(p % r for r in range(2, int(p ** .5) + 1))]
z = next((a, b) for a in range(q) for b in range(1, q) if all(pw((a, b), order // p) != (1, 0) for p in fac))
g = pw(z, q - 1)
T = [pw(g, k) for k in range(q + 1)]
m = q + 1
v = q + 2
def tr(x): return (2 * x[0]) % q
seqs = {}  # name -> (+-1 array, symmetric flag)
def add(name, vals):
    vals = np.array(vals)
    zeros = np.where(vals == 0)[0]
    for fills in itertools.product((1, -1), repeat=len(zeros)):
        w = vals.copy()
        w[zeros] = fills
        seqs[name + ('' if not len(zeros) else ' fill' + ''.join('+' if f > 0 else '-' for f in fills))] = w
delta = z  # element whose norm is a non-square (z primitive)
for c in range(q):
    add('chi(Tr t + %d)' % c, [chi(tr(t) + c) for t in T])
    add('chi(Tr(z t) + %d)' % c, [chi(tr(mul(delta, t)) + c) for t in T])
for lam in range(1, q):
    for bname, beta in (('lam', (lam, 0)), ('lam z', mul((lam, 0), delta))):
        add('eta(t - %s=%d)' % (bname, lam), [eta(((t[0] - beta[0]) % q, (t[1] - beta[1]) % q)) for t in T])
distinct = {}
for k, w in seqs.items():
    distinct.setdefault(tuple(w), k)
for w, k in list(distinct.items()):
    for j in range(1, m):  # shifts, for A
        distinct.setdefault(tuple(np.roll(w, j)), k + ' shift %d' % j)
names = [k for k in distinct.values()]
X = np.array(list(distinct.keys()))
sym = np.array([all(x[k] == x[(-k) % m] for k in range(m)) for x in X])
S = np.abs(np.fft.fft(X, axis=1)) ** 2
S = np.round(S[:, 1:], 6)
print('q', q, 'v', v, 'sequences', len(names), 'symmetric', sym.sum(), flush=True)
# rowsum condition: a^2 + 2 b^2 + d^2 = 4 with even row sums
rs = X.sum(1)
# hash D spectra
Dh = defaultdict(list)
for j in range(len(names)):
    Dh[S[j].tobytes()].append(j)
found = 0
for i in np.where(sym)[0]:
    for b in range(len(names)):
        if 'shift' in names[b]: continue
        need = np.round(4 * v - S[i] - 2 * S[b], 6)
        for d in Dh.get(need.tobytes(), []):
            if 'shift' in names[d]: continue
            if rs[i] ** 2 + 2 * rs[b] ** 2 + rs[d] ** 2 == 4:
                found += 1
                if found <= 10:
                    print('A =', names[i], '| B =', names[b], '| D =', names[d], '| row sums', rs[i], rs[b], rs[d])
print('triples (A symmetric, B and D up to shift)', found)

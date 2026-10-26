"""ModSpan (Howell echelon over Z/p^k) against brute force: span membership for ALL targets in (Z/q)^n for small q, n; and planted sums are never rejected."""
import sys, os, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
import modlattice as ML

def brute_span(rows, q, n):
    span = {tuple([0] * n)}
    frontier = [tuple([0] * n)]
    gens = [tuple(int(x) % q for x in r) for r in rows]
    while frontier:
        x = frontier.pop()
        for g in gens:
            y = tuple((a + b) % q for a, b in zip(x, g))
            if y not in span: span.add(y); frontier.append(y)
    return span

rng = np.random.default_rng(11)
ok = True; ntest = 0
for (p, k, n) in ((2, 3, 3), (2, 2, 4), (3, 2, 3), (2, 4, 2), (5, 1, 4), (3, 3, 2), (2, 1, 5)):
    q = p ** k
    for trial in range(40):
        m = int(rng.integers(1, 5))
        rows = rng.integers(0, q, size=(m, n)) * rng.choice([1, p, p * p if k >= 2 else 1], size=(m, 1))
        S = ML.ModSpan(p, k, n)
        for r in rows: S.insert(r)
        span = brute_span(rows, q, n)
        for t in itertools.product(range(q), repeat=n):
            a = S.contains(np.array(t)); b = t in span; ntest += 1
            if a != b: ok = False; print('MISMATCH', p, k, n, rows.tolist(), t, a, b)
        # bulk reduction: rows of random matrix reduce to zero iff in span
        X = rng.integers(0, q, size=(30, n)); R = S.reduce_rows(X)
        for x, r in zip(X, R):
            if (not np.any(r % q)) != (tuple(int(v) % q for v in x) in span): ok = False; print('REDUCE MISMATCH', p, k, n, rows.tolist(), x.tolist())
print('ModSpan vs brute force:', ntest, 'membership tests ->', 'OK' if ok else 'FAIL')
# planted: random integer vectors lists, a chosen sum must pass every modulus
bad = 0
for trial in range(30):
    n = 8
    Vs = [rng.integers(-30, 30, size=(int(rng.integers(5, 60)), n)) for _ in range(4)]
    pick = [int(rng.integers(0, len(v))) for v in Vs]
    T = sum(v[i] for v, i in zip(Vs, pick))
    res = ML.check(Vs, T, ((2, 3), (2, 5), (3, 3), (5, 2), (7, 1)))
    if not all(res.values()): bad += 1
print('planted sums rejected:', bad, '(must be 0)')
sys.exit(0 if ok and bad == 0 else 1)

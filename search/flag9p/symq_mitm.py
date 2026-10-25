"""Exhaustive search of L = [[X, W R], [W R, -conj(X)]] (order 2m, symmetric quaternary Hadamard iff
PAF_U + PAF_V + PAF_W = 2m delta, X = U + iV symmetric with complementary supports, W real circulant), by matching
PAF vectors (shifts 1..(m-1)/2) of X-part and W-part. Normalisation: W(0) = 1."""
import itertools, sys, numpy as np
m = int(sys.argv[1]); h = (m - 1) // 2
def paf(x): F = np.fft.fft(x); return np.rint(np.fft.ifft(np.abs(F) ** 2).real).astype(np.int64)[1:h + 1]
table = {}
for ch in itertools.product(range(4), repeat=h + 1):          # orbit 0 and orbits {t, -t}, t = 1..h: U+, U-, V+, V-
    U = np.zeros(m, np.int64); V = np.zeros(m, np.int64)
    for o, c in enumerate(ch):
        pos = [0] if o == 0 else [o, m - o]
        tgt = U if c < 2 else V; sgn = 1 if c % 2 == 0 else -1
        tgt[pos] = sgn
    key = tuple(paf(U) + paf(V))
    table.setdefault(key, []).append(ch)
hits = 0; ex = None
for wb in itertools.product((1, -1), repeat=m - 1):
    W = np.array((1,) + wb)
    key = tuple(-paf(W))
    if key in table:
        hits += len(table[key])
        if ex is None: ex = (table[key][0], W)
print(f"m={m} (order {2*m}): solutions {hits}", "example" if ex else "")
if ex:
    ch, W = ex
    U = np.zeros(m); V = np.zeros(m)
    for o, c in enumerate(ch):
        pos = [0] if o == 0 else [o, m - o]
        (U if c < 2 else V)[pos] = 1 if c % 2 == 0 else -1
    X = U + 1j * V
    circ = lambda x: np.array([np.roll(x, k) for k in range(m)]); R = np.eye(m)[::-1]
    L = np.block([[circ(X), circ(W) @ R], [circ(W) @ R, -circ(X).conj()]])
    print("  exact check: symmetric", np.allclose(L, L.T), " L conj(L) = 2m I", np.allclose(L @ L.conj(), 2 * m * np.eye(2 * m)),
          " entries in {±1, ±i}", np.allclose(np.abs(L), 1))

"""Small Goethals-Seidel and propus quadruples by meet in the middle (orders w <= 13), random picks."""
import numpy as np, itertools, random
def pafs(X):
    F = np.fft.fft(X, axis=1); return np.rint(np.fft.ifft(np.abs(F) ** 2, axis=1).real).astype(np.int64)
def all_seqs(w):
    return np.array([[1 - 2 * ((k >> i) & 1) for i in range(w)] for k in range(1 << w)], np.int64)
def gs_quads(w, count=5, seed=1):
    S = all_seqs(w); S = S[S[:, 0] == 1]; P = pafs(S)[:, 1:w // 2 + 1]
    rng = random.Random(seed); out = []
    idx = {}
    n = len(S)
    pairs_l = [(i, j) for i in range(n) for j in range(i, n)]
    rng.shuffle(pairs_l)
    for (i, j) in pairs_l[:200000]:
        idx.setdefault(tuple(P[i] + P[j]), []).append((i, j))
    for (i, j) in pairs_l[200000:2000000]:
        key = tuple(-(P[i] + P[j]))
        if key in idx:
            k, l = idx[key][0]
            out.append([S[k], S[l], S[i], S[j]])
            if len(out) >= count: break
    return out
def propus_quads(w, count=5, seed=1):
    S = all_seqs(w); S = S[S[:, 0] == 1]; P = pafs(S)[:, 1:w // 2 + 1]
    symm = [i for i in range(len(S)) if all(S[i][t] == S[i][(-t) % w] for t in range(w))]
    rng = random.Random(seed); out = []
    idx = {}
    for i in symm:
        for j in range(len(S)):
            idx.setdefault(tuple(P[i] + P[j]), []).append((i, j))      # A, D
    order = list(range(len(S))); rng.shuffle(order)
    for b in order:
        key = tuple(-2 * P[b])
        if key in idx:
            for (i, j) in idx[key][:3]:
                out.append([S[i], S[b], S[j]])
            if len(out) >= count: break
    return out

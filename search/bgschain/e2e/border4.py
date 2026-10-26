"""k = 4 border search for the spike-at-0 pattern on Q = Z2^2 (own implementation of the clique search of tiny_e2e.py): K (12x12), E (12x16), F (16x12)"""
import itertools, numpy as np
def gs_Q(ysq, k):
    Ym = [np.array([[int(y[g ^ h]) for h in range(k)] for g in range(k)]) for y in ysq]
    A_, B_, C_, D_ = Ym
    return np.block([[A_, B_, C_, D_], [-B_, A_, D_, -C_], [-C_, -D_, A_, B_], [-D_, C_, -B_, A_]])
def find_border(ysq, k=4, ntry=2000, seed=7):
    b = 4 * (k - 1)
    Mc = gs_Q(ysq, k)
    def balanced(m):
        out = []
        for c in itertools.combinations(range(m), m // 2):
            t = -np.ones(m, dtype=np.int64); t[list(c)] = 1; out.append(t)
        return np.array(out)
    Bal = balanced(k)
    ALL = np.array([np.concatenate([Bal[a_], Bal[b_], Bal[c_], Bal[d_]]) for a_, b_, c_, d_ in itertools.product(range(len(Bal)), repeat=4)])
    rng = np.random.default_rng(seed)
    def random_E():
        for _ in range(2000):
            cand = np.arange(len(ALL)); ch = []
            while len(ch) < b and len(cand):
                vv = ALL[cand[rng.integers(len(cand))]]; ch.append(vv); cand = cand[(ALL[cand] @ vv) == 0]
            if len(ch) == b: return np.array(ch)
    eps = np.array(list(itertools.product((1, -1), repeat=b)))
    for t in range(ntry):
        E = random_E()
        if E is None: continue
        Zm = Mc @ E.T
        ok = np.all(np.abs(eps @ Zm.T) == b, axis=1)
        S = eps[ok]
        if len(S) < b: continue
        G_ = (S @ S.T) == 0
        nS = len(S); adj = [sum(1 << j for j in range(nS) if G_[i, j]) for i in range(nS)]
        best = [None]
        def rec(cl, cand):
            if len(cl) == b: best[0] = list(cl); return True
            if bin(cand).count('1') < b - len(cl): return False
            while cand:
                i = cand.bit_length() - 1; cand &= ~(1 << i)
                if rec(cl + [i], cand & adj[i]): return True
            return False
        rec([], (1 << nS) - 1)
        if best[0] is not None:
            K = S[best[0]].T; F = -(Zm @ K) // b
            if np.all(np.abs(F) == 1):
                assert np.array_equal(E @ E.T, 4 * k * np.eye(b)) and np.array_equal(K @ K.T, b * np.eye(b)) and np.array_equal(Mc @ E.T + F @ K.T, np.zeros((4 * k, b)))
                return K, E, F
    return None

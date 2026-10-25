"""Circulant 4-splits for the x9 plug-in (split4.py).
GS: a Goethals-Seidel quadruple (A, B, C, D) of order w; each sequence's positions get a layer label 0..3 and each layer
(A_a, B_a, C_a, D_a) must satisfy sum PAF = w delta. Then the GS arrays of the layers are disjoint W(4w, w) summing to
the GS array of (A, B, C, D).
PROPUS: a propus quadruple (A, B, B, D) (A symmetric): labels on A symmetric, B used in both copies; each layer needs
PAF_A + 2 PAF_B + PAF_D = w delta; the propus arrays of the layers are disjoint symmetric W(4w, w)."""
import sys, time, numpy as np
from ortools.sat.python import cp_model
sys.path.insert(0, __import__('os').path.dirname(__file__))
from split4 import plug, is_weighing
def circ(x): return np.array([np.roll(x, k) for k in range(len(x))]).astype(np.int64)
def rev(X): return X[:, ::-1]
def gs_array(A, B, C, D):
    A, B, C, D = map(circ, (A, B, C, D))
    return np.block([[A, rev(B), rev(C), rev(D)], [-rev(B), A, rev(D.T), -rev(C.T)],
                     [-rev(C), -rev(D.T), A, rev(B.T)], [-rev(D), rev(C.T), -rev(B.T), A]])
def propus_array(A, B, D):
    A, B, D = map(circ, (A, B, D)); R = np.eye(A.shape[0], dtype=np.int64)[::-1]
    return np.block([[A, B, B, D], [B.T, D @ R, -A, -B @ R], [B.T, -A, -D @ R, B @ R], [D.T, -B @ R, B @ R, -A]])
def solve(kind, seqs, T=600, workers=1, fix=True, extra=None):
    w = len(seqs[0]); m = cp_model.CpModel()
    if kind == 'gs':
        names = [0, 1, 2, 3]; mult = {0: 1, 1: 1, 2: 1, 3: 1}
    else:
        names = [0, 1, 2]; mult = {0: 1, 1: 2, 2: 1}
    L = {}
    for s in names:
        for t in range(w):
            if kind == 'propus' and s == 0 and t > (w - t) % w:      # A symmetric: label(t) = label(-t)
                continue
            for a in range(4): L[s, t, a] = m.NewBoolVar('')
            m.AddExactlyOne([L[s, t, a] for a in range(4)])
    def lab(s, t, a):
        t %= w
        if kind == 'propus' and s == 0 and t > (w - t) % w: t = (w - t) % w
        return L[s, t, a]
    for a in range(4):
        m.Add(sum(mult[s] * lab(s, t, a) for s in names for t in range(w)) == w)
        for d in range(1, w // 2 + 1):
            terms = []
            for s in names:
                x = seqs[s]
                for t in range(w):
                    u = (t + d) % w
                    if kind == 'propus' and s == 0:
                        pass
                    y = m.NewBoolVar('')
                    m.AddBoolAnd([lab(s, t, a), lab(s, u, a)]).OnlyEnforceIf(y)
                    m.AddBoolOr([lab(s, t, a).Not(), lab(s, u, a).Not()]).OnlyEnforceIf(y.Not())
                    terms.append(int(mult[s] * x[t] * x[u]) * y)
            m.Add(sum(terms) == 0)
    if fix: m.Add(lab(0, 0, 0) == 1)
    s_ = cp_model.CpSolver(); s_.parameters.max_time_in_seconds = T; s_.parameters.num_workers = workers
    t0 = time.time(); r = s_.Solve(m)
    st = s_.StatusName(r)
    if r not in (cp_model.OPTIMAL, cp_model.FEASIBLE): return st, time.time() - t0, None
    labels = {s: [next(a for a in range(4) if s_.Value(lab(s, t, a))) for t in range(w)] for s in names}
    return st, time.time() - t0, labels
def verify(kind, seqs, labels):
    w = len(seqs[0])
    lay = lambda a: [np.where(np.array(labels[s]) == a, seqs[s], 0) for s in sorted(labels)]
    if kind == 'gs':
        K = gs_array(*seqs); parts = [gs_array(*lay(a)) for a in range(4)]
    else:
        K = propus_array(*seqs); parts = [propus_array(*lay(a)) for a in range(4)]
    N = K.shape[0]
    ok = np.array_equal(K @ K.T, N * np.eye(N, dtype=np.int64)) and np.array_equal(sum(parts), K) \
        and all(is_weighing(P, w) for P in parts) and all(np.array_equal(np.abs(P) * np.abs(Q), 0 * P) for i, P in enumerate(parts) for Q in parts[i + 1:])
    sym = kind == 'propus' and all((P == P.T).all() for P in parts)
    H = plug(parts); n = H.shape[0]
    had = np.array_equal(H @ H.T, n * np.eye(n, dtype=np.int64))
    return ok, sym, had, n

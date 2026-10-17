"""Exhaustive solution counts vs the Gaussian estimate on small Spence structures,
and an end-to-end Hadamard check of one solution per structure."""
import os
import itertools, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'verify'))
import verify
import spence_lib as L

def blockdiag(*blocks):
    n = sum(b.shape[0] for b in blocks)
    M = np.zeros((n, n), dtype=np.int64); i = 0
    for b in blocks:
        k = b.shape[0]; M[i:i+k, i:i+k] = b; i += k
    return M

def pool_vectors(st, r):
    out = []
    for bits in itertools.product((1.0, -1.0), repeat=st.O):
        y = np.array(bits)
        if int(y @ st.sizes) == r:
            out.append(y)
    return np.array(out).reshape(-1, st.O)

def exhaustive(st):
    v = st.G.v; pm = v % 4; target = (-4 - 4 * pm) // 4
    A, B = pool_vectors(st, 2), pool_vectors(st, 0)
    if len(A) == 0 or len(B) == 0:
        return 0, None, len(A), len(B)
    if len(A) * len(B) > 4e7 or len(B) ** 2 > 4e7:
        return -1, None, len(A), len(B)
    VA = np.rint((st.class_values(A) - pm) / 4).astype(np.int64)
    VB = np.rint((st.class_values(B) - pm) / 4).astype(np.int64)
    rng = np.random.default_rng(7)
    R = rng.integers(1, 2**62, size=st.C, dtype=np.int64).astype(np.uint64) | np.uint64(1)
    hA = (VA.astype(np.uint64) * R).sum(1, dtype=np.uint64)
    hB = (VB.astype(np.uint64) * R).sum(1, dtype=np.uint64)
    ht = (np.full(st.C, target, dtype=np.int64).astype(np.uint64) * R).sum(dtype=np.uint64)
    left = (hA[:, None] + hB[None, :]).ravel()
    right = (ht - (hB[:, None] + hB[None, :])).ravel()
    left.sort()
    lo = np.searchsorted(left, right, 'left'); hi = np.searchsorted(left, right, 'right')
    cnt = int((hi - lo).sum())
    first = None
    idx = np.flatnonzero(hi > lo)
    for q in idx[:50]:
        k, l = divmod(int(q), len(VB))
        need = target - VB[k] - VB[l]
        for i in range(len(VA)):
            m = np.flatnonzero((VB == (need - VA[i])).all(1))
            if len(m):
                first = (A[i], B[m[0]], B[k], B[l]); break
        if first is not None: break
    return cnt, first, len(A), len(B)

def to_cert(st, quad):
    seqs = []
    for y in quad:
        x = y[st.orb]
        seqs.append(''.join('+' if e > 0 else '-' for e in x))
    H4 = [[1,1,1,1],[1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]]
    P = (np.ones((4,4), int) - 2*np.eye(4, dtype=int)).tolist()
    return {"group": st.G.mods.tolist(), "sequences": seqs, "P": P, "Q": H4}

cases = []
I1 = np.eye(1, dtype=np.int64)
def I(k): return np.eye(k, dtype=np.int64)
S8, F8 = L.field_mult_matrix(2,3), L.frobenius_matrix(2,3)
S4 = L.field_mult_matrix(2,2)
S16 = L.field_mult_matrix(2,4)
S9 = L.field_mult_matrix(3,2)
def add(name, mods, mats): cases.append((name, L.AbGroup(mods), mats))
add("v=24  Z2^3xZ3   <F>", [2,2,2,3], [blockdiag(F8, I1)])
add("v=24  Z2^3xZ3   <F*(-1)>", [2,2,2,3], [blockdiag(F8, 2*I1)])
add("v=28  Z2^2xZ7   <S4*2>", [2,2,7], [blockdiag(S4, 2*I1)])
add("v=40  Z2^3xZ5   <S8>", [2,2,2,5], [blockdiag(S8, I1)])
add("v=40  Z2^3xZ5   <F>", [2,2,2,5], [blockdiag(F8, I1)])
add("v=48  Z2^4xZ3   <S16^5>", [2,2,2,2,3], [blockdiag(L.mat_pow_mod(S16,5,2), I1)])
add("v=48  Z2^4xZ3   <S16^3>", [2,2,2,2,3], [blockdiag(L.mat_pow_mod(S16,3,2), I1)])
add("v=56  Z2^3xZ7   <F*2>", [2,2,2,7], [blockdiag(F8, 2*I1)])
add("v=72  Z2^3xZ9   <S8,4>", [2,2,2,9], [blockdiag(S8, I1), blockdiag(I(3), 4*I1)])
add("v=72  Z2^3xZ3^2 <S8,S9^2>", [2,2,2,3,3], [blockdiag(S8, I(2)), blockdiag(I(3), L.mat_pow_mod(S9,2,3))])
add("v=88  Z2^3xZ11  <F,3>", [2,2,2,11], [blockdiag(F8, I1), blockdiag(I(3), 3*I1)])
add("v=104 Z2^3xZ13  <S8,3>", [2,2,2,13], [blockdiag(S8, I1), blockdiag(I(3), 3*I1)])
add("v=152 Z2^3xZ19  <S8,7>", [2,2,2,19], [blockdiag(S8, I1), blockdiag(I(3), 7*I1)])
add("v=104 Z2^3xZ13  <F,3>", [2,2,2,13], [blockdiag(F8, I1), blockdiag(I(3), 3*I1)])
add("v=120 Z2^3xZ15  <S8,4>", [2,2,2,15], [blockdiag(S8, I1), blockdiag(I(3), 4*I1)])
add("v=152 Z2^3xZ19  <F,7>", [2,2,2,19], [blockdiag(F8, I1), blockdiag(I(3), 7*I1)])
add("v=184 Z2^3xZ23  <S8,2>", [2,2,2,23], [blockdiag(S8, I1), blockdiag(I(3), 2*I1)])
for name, G, mats in cases:
    st = L.Structure(G, [G.perm_from_matrix(M) for M in mats], name)
    if st.O > 24:
        print(f"{name}: O={st.O} too many orbits for exhaustive"); continue
    est = st.estimate(samples=4000)
    cnt, first, nA, nB = exhaustive(st)
    lE = est['lE'] if est else float('-inf')
    msg = f"{name:28s} O={st.O:2d} C={st.C:2d} pools(2,0)=({nA},{nB}) exact={cnt:>10d}  Gaussian=10^{lE:.2f}"
    if first is not None:
        H, _ = verify.build_spence(to_cert(st, first))
        n = H.shape[0]
        ok = np.isin(H, (-1,1)).all() and np.array_equal(H @ H.T, n*np.eye(n, dtype=int))
        msg += f"  -> H({n}) Hadamard: {ok}"
    print(msg, flush=True)

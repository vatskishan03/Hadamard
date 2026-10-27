import numpy as np, itertools, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import cmb, gsplug as GP
U4 = np.array([1, 1j, -1, -1j])
def negacyclic(x):
    m = len(x); M = np.zeros((m, m), complex)
    for i in range(m):
        for j in range(m):
            M[i, j] = x[j - i] if j >= i else -x[m + j - i]
    return M
def find_pair(m, count=1):
    mats = {}
    for e in itertools.product(range(4), repeat=m):
        x = U4[list(e)]
        M = negacyclic(x)
        mats[e] = M
    # AA* + BB* = 2m I : compare first rows of A A* (negacyclic, hermitian-ish)
    keyd = {}
    for e, M in mats.items():
        G = M @ M.conj().T
        keyd.setdefault(tuple(np.round(G[0], 6)), []).append(e)
    out = []
    target = np.zeros(m, complex); target[0] = 2 * m
    for k, es in keyd.items():
        k2 = tuple(np.round(target - np.array(k), 6))
        if k2 in keyd:
            out.append((es[0], keyd[k2][0]))
            if len(out) >= count: return out
    return out
def reversal_nega(q):
    m = (q - 1) // 2
    R = np.arange(q)
    for h in range(2):
        for j in range(m):
            R[1 + h * m + j] = 1 + (h if j == 0 else 1 - h) * m + (-j) % m
    return R
for m in (3, 5):
    q = 2 * m + 1
    pr = find_pair(m)
    if not pr:
        print("m =", m, "no negacyclic pair"); continue
    a, b = pr[0]
    A, B = negacyclic(U4[list(a)]), negacyclic(U4[list(b)])
    K = np.block([[A, B], [-B.conj().T, A.conj().T]])
    assert np.allclose(K @ K.conj().T, 2 * m * np.eye(2 * m))
    Kr, Ki = np.rint(K.real).astype(np.int8), np.rint(K.imag).astype(np.int8)
    N = cmb.build_blocks(q, cmb.prim_root(q), Kr, Ki)
    # patch reversal in check_blocks
    cmb.reversal_perm = lambda qq: reversal_nega(qq)
    GP.cmb.reversal_perm = cmb.reversal_perm
    print("m", m, "negacyclic pair", a, b, ":", end=" ")
    ok = cmb.check_blocks(N)
    Bk = GP.real_blocks(N)
    print("  real blocks:", GP.check_real_blocks(Bk, q))
    for t in (3, 5):
        T = cmb.t_matrices(t)
        Aj = GP.plug_A(Bk, T)
        Pm, Sm = GP.gs_R(q, t)
        G = GP.gs_array(Aj, Pm, Sm)
        print(f"   q={q} t={t}: order {G.shape[0]}, H H^T = nI:", GP.is_hadamard(G))

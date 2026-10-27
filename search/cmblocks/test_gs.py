import numpy as np, sys, time
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import cmb, gsplug as GP, paley_circ as P
for qP in (11, 19, 43):
    Kr, Ki, m = P.block_circulant_K(qP)
    q = 2 * m + 1
    g = cmb.prim_root(q)
    N = cmb.build_blocks(q, g, Kr, Ki)
    Bk = GP.real_blocks(N)
    print("q", q, "real blocks: R'BR'=B^T, commute, orthogonal rows:", GP.check_real_blocks(Bk, q))
    for t in (3, 5):
        T = cmb.t_matrices(t)
        A = GP.plug_A(Bk, T)
        Pm, Sm = GP.gs_R(q, t)
        assert all(np.abs(x).max() == 1 and np.abs(x).min() == 1 for x in A)
        G = GP.gs_array(A, Pm, Sm)
        n = G.shape[0]
        print(f"  q={q} t={t}: order {n} = 8qt, entries +-1: {np.all(np.abs(G)==1)}, H H^T = nI: {GP.is_hadamard(G)}")

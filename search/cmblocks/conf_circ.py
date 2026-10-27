"""Paley's symmetric conference matrix C of order p+1 (p = 1 mod 4 prime) and K = C + iI (quaternary Hadamard of order
p+1 = 2N, N odd) in exact 2 x 2 block circulant form: rho = tau^2 (order N = (p+1)/2) in the non-split torus acts
semiregularly on PG(1,p) with two orbits (monomial action P_rho commutes with C; gauge E_{a,j} = P_rho^j e_{f_a})."""
import numpy as np, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import paley_circ as P

def block_circulant_K_conf(p, check=True):
    assert p % 4 == 1 and P.is_prime(p)
    n = p + 1; N = n // 2
    chi = P.chi_table(p)
    C = P.paley_H(p, chi) - np.eye(n, dtype=np.int64)
    assert np.array_equal(C, C.T) and np.array_equal(C @ C.T, p * np.eye(n, dtype=np.int64))
    tau = P.find_tau(p)
    it, st = P.action(tau, p, chi)
    def power(k):
        img, sg = np.arange(n), np.ones(n, np.int64)
        bi, bs = it.copy(), st.copy()
        while k:
            if k & 1: img, sg = P.compose(bi, bs, img, sg)
            bi, bs = P.compose(bi, bs, bi, bs)
            k >>= 1
        return img, sg
    ri, rs = power(2)                       # rho = tau^2, order N
    ci, cs = power(N * 2)
    assert np.array_equal(ci, np.arange(n)) and np.all(cs == 1)
    O0 = set(); u = 0
    while u not in O0:
        O0.add(u); u = int(ri[u])
    assert len(O0) == N
    f1 = next(u for u in range(n) if u not in O0)
    pts = np.empty(2 * N, np.int64); sig = np.empty(2 * N, np.int64)
    for a, f in enumerate((0, f1)):
        u, s = f, 1
        for j in range(N):
            pts[a * N + j], sig[a * N + j] = u, s
            s *= rs[u]; u = int(ri[u])
        assert u == f and s == 1
    Kr = sig[:, None] * sig[None, :] * C[np.ix_(pts, pts)]
    Ki = np.eye(2 * N, dtype=np.int64)
    if check:
        Kc = Kr + 1j * Ki
        assert np.allclose(Kc @ Kc.conj().T, n * np.eye(n))
        for A in range(2):
            for B in range(2):
                blk = Kc[A * N:(A + 1) * N, B * N:(B + 1) * N]
                for j in range(N):
                    assert np.array_equal(blk[j], np.roll(blk[0], j))
    return Kr.astype(np.int8), Ki.astype(np.int8), N

if __name__ == "__main__":
    import cmb, gsplug as GP
    for p in (5, 17, 29, 41):
        Kr, Ki, N = block_circulant_K_conf(p)
        q = 2 * N + 1
        assert q == p + 2 and P.is_prime(q)
        Nb = cmb.build_blocks(q, cmb.prim_root(q), Kr, Ki)
        ok = cmb.check_blocks(Nb, verbose=False)
        Bk = GP.real_blocks(Nb)
        res = GP.check_real_blocks(Bk, q)
        out = []
        for t in (3,):
            G = GP.gs_array(GP.plug_A(Bk, cmb.t_matrices(t)), *GP.gs_R(q, t))
            out.append((t, G.shape[0], GP.is_hadamard(G)))
        print(f"conference p={p}: K = C + iI of order {2*N} block circulant; Miyamoto q={q}: blocks ok {ok}, real {res}, GS {out}")

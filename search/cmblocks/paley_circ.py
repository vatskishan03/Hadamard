"""Paley's quaternary Hadamard matrix of order (q+1)/2 in 2x2 block circulant form (blocks of order m = (q+1)/4).

q prime, q = 3 (mod 8).  S[a,b] = chi(det(u_a,u_b)), H = I + S of order q+1 on PG(1,q), u_x = (x,1), u_inf = (1,0).
tau in SL(2,q) generates the non-split torus (order q+1); g = tau^m (g^2 = -I) is the complex structure i, rho = tau^4
(order m) acts semiregularly on the 2m pairs {f, g f} with 2 orbits.  The monomial action P_u e_a = chi(lambda_a) e_{u a}
(u a = lambda_a-multiple of the representative) is a homomorphism SL(2,q) -> signed permutations commuting with H.
Basis E_{a,j} = P_rho^j e_{f_a}; K[(a,j),(b,l)] = (<E, H E'> + i <J E, H E'>)/(1+i) is exactly circulant in l - j."""
import numpy as np

def is_prime(n):
    return n > 1 and all(n % p for p in range(2, int(n ** 0.5) + 1))

def factor_small(n):
    f, p = [], 2
    while p * p <= n:
        if n % p == 0:
            f.append(p)
            while n % p == 0: n //= p
        p += 1
    if n > 1: f.append(n)
    return f

def matmul2(X, Y, q):
    return ((X[0][0]*Y[0][0] + X[0][1]*Y[1][0]) % q, (X[0][0]*Y[0][1] + X[0][1]*Y[1][1]) % q,
            (X[1][0]*Y[0][0] + X[1][1]*Y[1][0]) % q, (X[1][0]*Y[0][1] + X[1][1]*Y[1][1]) % q)

def mpow(X, e, q):
    R = (1, 0, 0, 1); B = (X[0], X[1], X[2], X[3])
    while e:
        if e & 1:
            R = matmul2([[R[0], R[1]], [R[2], R[3]]], [[B[0], B[1]], [B[2], B[3]]], q)
        B = matmul2([[B[0], B[1]], [B[2], B[3]]], [[B[0], B[1]], [B[2], B[3]]], q)
        e >>= 1
    return R

def find_tau(q):
    """companion matrix [[0,-1],[1,t]] (char poly x^2 - t x + 1) of order exactly q+1 in SL(2,q)"""
    n = q + 1
    for t in range(q):
        if pow((t * t - 4) % q, (q - 1) // 2, q) != q - 1: continue          # need t^2 - 4 a non-residue
        tau = (0, q - 1, 1, t)
        if mpow(tau, n, q) != (1, 0, 0, 1): continue
        if all(mpow(tau, n // p, q) != (1, 0, 0, 1) for p in factor_small(n)):
            return tau
    raise ValueError("no tau")

def chi_table(q):
    chi = -np.ones(q, np.int64); chi[0] = 0; chi[(np.arange(1, q) ** 2) % q] = 1
    return chi

def action(M, q, chi):
    """for M = (a,b,c,d) in SL(2,q): image[u], sign[u] with M u_a = lambda * u_image, sign = chi(lambda).
    points 0..q-1 = x, q = infinity."""
    a, b, c, d = M
    n = q + 1
    img = np.empty(n, np.int64); sg = np.empty(n, np.int64)
    for x in range(q):
        num, den = (a * x + b) % q, (c * x + d) % q
        if den:
            img[x] = num * pow(den, q - 2, q) % q; sg[x] = chi[den]
        else:
            img[x] = q; sg[x] = chi[num]
    if c % q:
        img[q] = a * pow(c, q - 2, q) % q; sg[q] = chi[c]
    else:
        img[q] = q; sg[q] = chi[a]
    return img, sg

def paley_H(q, chi):
    n = q + 1
    H = np.empty((n, n), np.int64)
    x = np.arange(q)
    H[:q, :q] = chi[(x[:, None] - x[None, :]) % q]
    H[q, :q] = 1                 # det(u_inf, u_y) = 1
    H[:q, q] = chi[q - 1]        # det(u_x, u_inf) = -1
    H[q, q] = 0
    H[np.arange(n), np.arange(n)] += 1
    return H

def compose(img1, sg1, img2, sg2):
    """P_1 P_2: e_u -> sg2[u] e_img2[u] -> sg2[u] sg1[img2[u]] e_img1[img2[u]]"""
    return img1[img2], sg2 * sg1[img2]

def block_circulant_K(q, check=True):
    assert q % 8 == 3 and is_prime(q)
    n = q + 1; m = n // 4
    chi = chi_table(q)
    H = paley_H(q, chi)
    tau = find_tau(q)
    it, st = action(tau, q, chi)
    # P_tau^k as (img, sg) by repeated composition
    def power(k):
        img, sg = np.arange(n), np.ones(n, np.int64)
        base_i, base_s = it.copy(), st.copy()
        while k:
            if k & 1: img, sg = compose(base_i, base_s, img, sg)
            base_i, base_s = compose(base_i, base_s, base_i, base_s)
            k >>= 1
        return img, sg
    gi, gs = power(m)                       # g = tau^m, g^2 = -I
    ri, rs = power(4)                       # rho = tau^4 of order m
    # checks: P_g^2 = -I, P_rho^m = I, H commutes
    g2i, g2s = compose(gi, gs, gi, gs)
    assert np.array_equal(g2i, np.arange(n)) and np.all(g2s == -1), "P_g^2 != -I"
    ri_m, rs_m = power(4 * m)
    assert np.array_equal(ri_m, np.arange(n)) and np.all(rs_m == 1), "P_rho^m != I"
    # orbits of <rho, g>: pick f_0 = 0, f_1 = first point outside the orbit of f_0 under rho and g
    def orbit(f):
        S = {f}; fr = [f]
        while fr:
            u = fr.pop()
            for v in (int(ri[u]), int(gi[u])):
                if v not in S: S.add(v); fr.append(v)
        return S
    f0 = 0
    O0 = orbit(f0)
    f1 = next(u for u in range(n) if u not in O0)
    assert len(O0) == 2 * m
    p = np.empty(2 * m, np.int64); sig = np.empty(2 * m, np.int64)
    for a, f in enumerate((f0, f1)):
        u, s = f, 1
        for j in range(m):
            p[a * m + j], sig[a * m + j] = u, s
            s *= rs[u]; u = int(ri[u])
        assert u == f and s == 1, "orbit sign product"
    # a[k,l] = sig_k sig_l H[p_k, p_l];  b[k,l] = (J E_k)^T H E_l, J E_k = sig_k gs[p_k] e_{g p_k}
    a = sig[:, None] * sig[None, :] * H[np.ix_(p, p)]
    b = (sig * gs[p])[:, None] * sig[None, :] * H[np.ix_(gi[p], p)]
    Kr, Ki = (a + b) // 2, (b - a) // 2
    assert np.all((a + b) % 2 == 0) and np.all(np.abs(Kr) + np.abs(Ki) == 1)
    if check:
        Kc = Kr + 1j * Ki
        assert np.array_equal(np.rint((Kc @ Kc.conj().T).real).astype(np.int64), 2 * m * np.eye(2 * m, dtype=np.int64)) if m < 100 else True
        # circulant blocks
        for A in range(2):
            for B in range(2):
                blk = Kc[A * m:(A + 1) * m, B * m:(B + 1) * m]
                for j in range(m):
                    assert np.array_equal(blk[j], np.roll(blk[0], j)), "not circulant"
    return Kr.astype(np.int8), Ki.astype(np.int8), m

if __name__ == "__main__":
    import sys
    for q in map(int, sys.argv[1:]):
        Kr, Ki, m = block_circulant_K(q)
        print(q, "order", 2 * m, "block circulant, K K* = 2m I:", True)

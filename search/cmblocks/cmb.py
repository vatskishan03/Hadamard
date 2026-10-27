"""Quaternary Miyamoto blocks (q = 3 mod 4 prime) and the real Goethals-Seidel plug-in on the realification.

Blocks:  N_rc = [[A_rc, B_rc e], [B_rc e^T, [[U+V, U-V], [U-V, U+V]]_rc]]  (order q = 1 + 2m), borders (A, B) a 4x4 pair
found by border_search.py with  A - 2R = B  (R the block row sums of the Hermitian core), B B* = A A* = 4I,
A B* + 2 B R* = 0.  The core W = i [chi(y-x)] is listed along a primitive root g (squares g^(2j), then non-squares
g^(2j+1)), K = 2x2 block circulant quaternary Hadamard matrix of order 2m = q - 1, V = [[I, K], [-K*, I]]."""
import itertools, numpy as np, sys

# ---------------- T-matrices (Turyn: base sequences) ----------------
def base_sequences(m):
    def nac(x):
        return tuple(int(np.dot(x[:-k], x[k:])) for k in range(1, len(x)))
    seqs = lambda L: [np.array(b) for b in itertools.product((1, -1), repeat=L)]
    S1, S0 = seqs(m + 1), seqs(m)
    need = {}
    for c in S0:
        for d in S0:
            v = tuple(a + b for a, b in zip(nac(c), nac(d))) + (0,) if m > 0 else (0,)
            need.setdefault(v, (c, d))
    for a in S1:
        for b in S1:
            v = tuple(-(x + y) for x, y in zip(nac(a), nac(b)))
            key = v[:-1] + (0,) if m > 0 else (0,)
            if v[-1] == 0 and key in need:
                c, d = need[key]
                return a, b, c, d
    raise ValueError("no base sequences")

def t_matrices(t):
    m = (t - 1) // 2
    if m == 0:
        return np.array([[1], [0], [0], [0]], np.int64)
    a, b, c, d = base_sequences(m)
    z1, z0 = np.zeros(m + 1, np.int64), np.zeros(m, np.int64)
    T = np.array([np.concatenate([(a + b) // 2, z0]), np.concatenate([(a - b) // 2, z0]),
                  np.concatenate([z1, (c + d) // 2]), np.concatenate([z1, (c - d) // 2])], np.int64)
    check_t(T)
    return T

def circ_int(v):
    v = np.asarray(v); t = len(v)
    idx = (np.arange(t)[None, :] - np.arange(t)[:, None]) % t
    return v[idx]                      # row i = v rolled by i (row 0 = v)

def check_t(T):
    T = np.asarray(T); t = T.shape[1]
    assert T.shape[0] == 4 and (np.abs(T).sum(0) == 1).all()
    assert np.array_equal(sum(circ_int(r) @ circ_int(r).T for r in T), t * np.eye(t, dtype=np.int64))

# ---------------- primitive root, characters ----------------
def factor_set(n):
    f, p = set(), 2
    while p * p <= n:
        while n % p == 0: f.add(p); n //= p
        p += 1
    if n > 1: f.add(n)
    return f

def prim_root(q):
    ps = factor_set(q - 1)
    return next(g for g in range(2, q) if all(pow(g, (q - 1) // p, q) != 1 for p in ps))

def chi_table(q):
    chi = -np.ones(q, np.int64); chi[0] = 0; chi[(np.arange(1, q) ** 2) % q] = 1
    return chi

# ---------------- borders ----------------
def border_pair():
    """(A, B), first solution of border_search.py for R = [[0,i,0,0],[-i,0,0,0],[0,0,0,i],[0,0,-i,0]], lam = kappa = 1"""
    sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))
    import border_search as bs
    R = bs.variants()[(1, 1)]
    lam, kappa, B, A = bs.search(R)[0]
    assert lam == 1 and kappa == 1
    return np.rint(A.real).astype(np.int64) + 1j * np.rint(A.imag).astype(np.int64), \
           np.rint(B.real).astype(np.int64) + 1j * np.rint(B.imag).astype(np.int64)

def core_blocks(q, g, Kr, Ki):
    """U (4m x 4m, complex) and V (4m x 4m) as 4x4 arrays of m x m blocks"""
    m = (q - 1) // 2
    chi = chi_table(q)
    pts = np.array([pow(g, 2 * j, q) for j in range(m)] + [pow(g, 2 * j + 1, q) for j in range(m)], np.int64)
    Wi = chi[(pts[None, :] - pts[:, None]) % q]                    # W = i Wi, Wi[x, y] = chi(y - x)
    K = Kr.astype(np.int64) + 1j * Ki.astype(np.int64)
    I2 = np.eye(2 * m)
    V = np.block([[I2, K], [-K.conj().T, I2]])
    U = 1j * np.block([[Wi, np.zeros_like(Wi)], [np.zeros_like(Wi), Wi]])
    blk = lambda X, r, c: X[r * m:(r + 1) * m, c * m:(c + 1) * m]
    return [[blk(U, r, c) for c in range(4)] for r in range(4)], [[blk(V, r, c) for c in range(4)] for r in range(4)]

def build_blocks(q, g, Kr, Ki, A=None, B=None):
    m = (q - 1) // 2
    if A is None: A, B = border_pair()
    U, V = core_blocks(q, g, Kr, Ki)
    N = [[None] * 4 for _ in range(4)]
    e = np.ones(2 * m)
    for r in range(4):
        for c in range(4):
            Nrc = np.zeros((q, q), complex)
            Nrc[0, 0] = A[r, c]
            Nrc[0, 1:] = B[r, c]
            Nrc[1:, 0] = B[r, c]
            u, v = U[r][c], V[r][c]
            Nrc[1:m + 1, 1:m + 1] = u + v; Nrc[1 + m:, 1 + m:] = u + v
            Nrc[1:m + 1, 1 + m:] = u - v;  Nrc[1 + m:, 1:m + 1] = u - v
            N[r][c] = Nrc
    return N

def assemble(N):
    q = N[0][0].shape[0]
    return np.block(N)

def reversal_perm(q):
    """R_B: fixes the border, (h, j) -> (h, -j)"""
    m = (q - 1) // 2
    R = np.arange(q)
    for h in range(2):
        for j in range(m):
            R[1 + h * m + j] = 1 + h * m + (-j) % m
    return R

def check_blocks(N, verbose=True):
    q = N[0][0].shape[0]
    H = assemble(N)
    n = H.shape[0]
    ok_unit = np.all(np.isclose(np.abs(H), 1)) and np.all(np.isclose(H.real, np.rint(H.real))) and np.all(np.isclose(H.imag, np.rint(H.imag)))
    G = H @ H.conj().T
    ok_had = np.allclose(G, n * np.eye(n))
    Rb = reversal_perm(q)
    flat = [N[r][c] for r in range(4) for c in range(4)]
    ok_rt = all(np.array_equal(x[np.ix_(Rb, Rb)], x.T) for x in flat)
    ok_comm = all(np.allclose(x @ y, y @ x) for i, x in enumerate(flat) for y in flat[i:])
    if verbose:
        print(f"q={q}: entries units {ok_unit}, H H* = 4qI {ok_had}, R N R = N^T {ok_rt}, blocks commute {ok_comm}")
    return ok_unit and ok_had and ok_rt and ok_comm

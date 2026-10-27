#!/usr/bin/env python3
"""Rebuild and check the Hadamard matrix of order 8*q*t of a 'quaternary Miyamoto blocks x T-matrices' certificate.

    python3 rebuild.py CERT.json [--no-gram]      (needs only numpy; memory about 1 GB for q = 863, t = 3)

Ingredients of the certificate (all exact integer data):
  q          prime = 3 (mod 4), m = (q - 1)/2 odd, g a primitive root mod q
  K          quaternary Hadamard matrix of order 2m = q - 1, as 2 x 2 blocks of order m; each block is circulant, and
             the certificate gives the first row of each block as digits e (entry i^e, i = sqrt(-1))
  A, B       the 4 x 4 border matrices of units (digits e); R = rows sums of the core, see MIYAMOTO_BLOCKS_Q.md
  T          four first rows over + - 0 of circulant T-matrices of order t
Construction (all steps checked here):
  1. W = i [chi(y - x)] on the points g^(2j), g^(2j+1) (j in Z_m); U = diag(W, W); V = [[I, K], [-K*, I]] (order 4m);
     sixteen blocks of order q, N_rc = [[A_rc, B_rc e], [B_rc e^T, [[U+V, U-V], [U-V, U+V]]_rc]]  (U_rc, V_rc of order m);
     H = (N_rc) is a quaternary Hadamard matrix of order 4q (H H* = 4q I).
  2. X_rc = Re N_rc + Im N_rc, Y_rc = Re N_rc - Im N_rc;  B_rc = [[X_rc, Y_rc], [-Y_rc, X_rc]] (order 2q, +-1 entries):
     the sixteen real blocks commute and R' B_rc R' = B_rc^T for the signed permutation R' = diag(1, -1) (x) R_B,
     R_B: (h, j) -> (h, -j), border fixed.
  3. A_j = sum_k T_k (x) B_kj (order 2qt), R'' = R_t (x) R' (R_t: s -> -s on Z_t) and the Goethals-Seidel array
        [[A1, A2R, A3R, A4R], [-A2R, A1, A4^T R, -A3^T R], [-A3R, -A4^T R, A1, A2^T R], [-A4R, A3^T R, -A2^T R, A1]]
     is a Hadamard matrix of order 8qt: G G^T = 8qt I.
"""
import hashlib, json, sys
import numpy as np

U4 = np.array([1, 1j, -1, -1j])

def circulant_rows(first):
    m = len(first)
    return first[(np.arange(m)[None, :] - np.arange(m)[:, None]) % m]

def build_H(cert):
    q, g, m = cert["q"], cert["g"], (cert["q"] - 1) // 2
    assert q % 4 == 3 and all(q % p for p in range(2, int(q ** 0.5) + 1))
    assert len({pow(g, k, q) for k in range(q - 1)}) == q - 1, "g is not a primitive root"
    rows = cert["K"]["rows"] if "rows" in cert["K"] else cert["K"]
    Kb = [[U4[np.array([int(c) for c in rows[f"{a}{b}"]])] for b in range(2)] for a in range(2)]
    K = np.block([[circulant_rows(Kb[a][b]) for b in range(2)] for a in range(2)])
    assert K.shape == (2 * m, 2 * m)
    chi = -np.ones(q, np.int64); chi[0] = 0; chi[(np.arange(1, q) ** 2) % q] = 1
    pts = np.array([pow(g, 2 * j, q) for j in range(m)] + [pow(g, 2 * j + 1, q) for j in range(m)], np.int64)
    Wi = chi[(pts[None, :] - pts[:, None]) % q]
    Z = np.zeros_like(Wi)
    U = 1j * np.block([[Wi, Z], [Z, Wi]])
    I2 = np.eye(2 * m)
    V = np.block([[I2, K], [-K.conj().T, I2]])
    A = U4[np.array(cert["A"])]; B = U4[np.array(cert["B"])]
    blk = lambda X, r, c: X[r * m:(r + 1) * m, c * m:(c + 1) * m]
    N = [[None] * 4 for _ in range(4)]
    for r in range(4):
        for c in range(4):
            Nrc = np.zeros((q, q), complex)
            Nrc[0, 0] = A[r, c]; Nrc[0, 1:] = B[r, c]; Nrc[1:, 0] = B[r, c]
            u, v = blk(U, r, c), blk(V, r, c)
            Nrc[1:m + 1, 1:m + 1] = u + v; Nrc[1 + m:, 1 + m:] = u + v
            Nrc[1:m + 1, 1 + m:] = u - v; Nrc[1 + m:, 1:m + 1] = u - v
            N[r][c] = Nrc
    return N

def quaternary_check(N):
    """H H* = 4q I for H = (N_rc), exactly (float32 products of small integers)"""
    q = N[0][0].shape[0]
    H = np.block(N)
    Hr, Hi = np.rint(H.real).astype(np.float32), np.rint(H.imag).astype(np.float32)
    assert np.all(np.abs(Hr) + np.abs(Hi) == 1), "entries are not units"
    n = H.shape[0]
    re = Hr @ Hr.T + Hi @ Hi.T
    im = Hi @ Hr.T - Hr @ Hi.T
    return np.array_equal(re, 4 * q * np.eye(n, dtype=np.float32)) and not im.any()

def real_blocks(N):
    out = [[None] * 4 for _ in range(4)]
    for r in range(4):
        for c in range(4):
            Hr, Hi = np.rint(N[r][c].real).astype(np.int8), np.rint(N[r][c].imag).astype(np.int8)
            X, Y = Hr + Hi, Hr - Hi
            out[r][c] = np.block([[X, Y], [-Y, X]]).astype(np.int8)
    return out

def reversal(q):
    m = (q - 1) // 2
    R = np.arange(q)
    for h in range(2):
        for j in range(m):
            R[1 + h * m + j] = 1 + h * m + (-j) % m
    return R

def Rprime(q):
    Rb = reversal(q)
    return np.concatenate([Rb, q + Rb]), np.concatenate([np.ones(q, np.int8), -np.ones(q, np.int8)])

def algebra_check(Bk, q):
    """sixteen real blocks of order 2q: R' B R' = B^T and pairwise commuting (float32, exact)"""
    perm, sign = Rprime(q)
    n = 2 * q
    ok_rt = True
    for r in range(4):
        for c in range(4):
            X = Bk[r][c]
            RXR = np.empty_like(X)
            RXR[np.ix_(perm, perm)] = X * sign[:, None] * sign[None, :]      # (R' X R')[perm[i], perm[j]] = s_i s_j X[i, j]
            ok_rt &= np.array_equal(RXR, X.T)
    F = [Bk[r][c].astype(np.float32) for r in range(4) for c in range(4)]
    ok_comm = True
    for i in range(16):
        for j in range(i + 1, 16):
            ok_comm &= np.array_equal(F[i] @ F[j], F[j] @ F[i])
    return ok_rt, ok_comm

def plug(Bk, T, q):
    t = T.shape[1]
    Ct = [circulant_rows(np.asarray(r)).astype(np.int8) for r in T]
    A = [sum(np.kron(Ct[k], Bk[k][j]) for k in range(4)).astype(np.int8) for j in range(4)]
    assert all(np.all(np.abs(x) == 1) for x in A)
    perm, sign = Rprime(q)
    n2 = 2 * q
    P = np.empty(t * n2, np.int64); S = np.empty(t * n2, np.int8)
    for s in range(t):
        P[s * n2:(s + 1) * n2] = ((-s) % t) * n2 + perm
        S[s * n2:(s + 1) * n2] = sign
    return A, P, S

def mulR(X, P, S):
    Y = np.empty_like(X)
    Y[:, P] = X * S[None, :]
    return Y

def gs_array(A, P, S):
    n = A[0].shape[0]
    AR = [mulR(x, P, S) for x in A]
    ATR = [mulR(np.ascontiguousarray(x.T), P, S) for x in A]
    rows = [[A[0], AR[1], AR[2], AR[3]], [-AR[1], A[0], ATR[3], -ATR[2]],
            [-AR[2], -ATR[3], A[0], ATR[1]], [-AR[3], ATR[2], -ATR[1], A[0]]]
    G = np.empty((4 * n, 4 * n), np.int8)
    for i in range(4):
        for j in range(4):
            G[i * n:(i + 1) * n, j * n:(j + 1) * n] = rows[i][j]
    return G

def gram_scalar(G, c, block=2048):
    """G G^T = c I exactly: float32 products of +-1 data, sums below 2^24, blocks on and above the diagonal"""
    n = G.shape[0]
    assert n < 2 ** 24
    for i in range(0, n, block):
        Fi = G[i:i + block].astype(np.float32)
        for j in range(i, n, block):
            M = Fi @ (Fi if j == i else G[j:j + block].astype(np.float32)).T
            if j == i:
                M[np.diag_indices_from(M)] -= c
            if M.any():
                return False
    return True

def csv_sha256(G, rows=512):
    h = hashlib.sha256()
    for i in range(0, G.shape[0], rows):
        X = G[i:i + rows]
        cell = np.empty(X.shape + (3,), np.uint8)
        cell[..., 0] = np.where(X < 0, ord("-"), 0); cell[..., 1] = ord("1"); cell[..., 2] = ord(",")
        cell[:, -1, 2] = ord("\n")
        cell = cell.ravel()
        h.update(cell[cell != 0].tobytes())
    return h.hexdigest()

def csv_stream(G, out, rows=512):
    """write the canonical CSV (rows of 1/-1 joined by commas, one per line) to the binary stream out"""
    for i in range(0, G.shape[0], rows):
        X = G[i:i + rows]
        cell = np.empty(X.shape + (3,), np.uint8)
        cell[..., 0] = np.where(X < 0, ord("-"), 0); cell[..., 1] = ord("1"); cell[..., 2] = ord(",")
        cell[:, -1, 2] = ord("\n")
        cell = cell.ravel()
        out.write(cell[cell != 0].tobytes())

def main(path, gram=True, save=None):
    cert = json.load(open(path))
    q = cert["q"]
    t = len(cert["T"][0])
    T = np.array([[{"+": 1, "-": -1, "0": 0}[c] for c in r] for r in cert["T"]], np.int64)
    Ct = [circulant_rows(r) for r in T]
    assert (np.abs(T).sum(0) == 1).all() and np.array_equal(sum(c @ c.T for c in Ct), t * np.eye(t, dtype=np.int64)), "not T-matrices"
    print(f"q = {q}, t = {t}, order 8qt = {8 * q * t}")
    N = build_H(cert)
    print("quaternary H of order 4q: H H* = 4q I :", quaternary_check(N))
    Bk = real_blocks(N)
    del N
    print("real blocks of order 2q: R' B R' = B^T, pairwise commuting :", algebra_check(Bk, q))
    A, P, S = plug(Bk, T, q)
    G = gs_array(A, P, S)
    n = G.shape[0]
    assert n == 8 * q * t
    print("Goethals-Seidel array: order", n, "entries +-1:", bool(np.all(np.abs(G) == 1)))
    if gram:
        print("G G^T = nI exactly:", gram_scalar(G, n))
    sha = csv_sha256(G)
    print("sha256 of the CSV rows:", sha)
    if save:
        np.save(save, G)
        print("saved", save)
    if "sha256" in cert:
        print("matches certificate:", sha == cert["sha256"])
    return sha

if __name__ == "__main__":
    save = sys.argv[sys.argv.index("--save") + 1] if "--save" in sys.argv else None
    if "--csv" in sys.argv:      # stream the CSV of the matrix to stdout (e.g. into check_matrix); no other output
        cert = json.load(open(sys.argv[1]))
        q = cert["q"]
        T = np.array([[{"+": 1, "-": -1, "0": 0}[c] for c in r] for r in cert["T"]], np.int64)
        Bk = real_blocks(build_H(cert))
        A, P, S = plug(Bk, T, q)
        csv_stream(gs_array(A, P, S), sys.stdout.buffer)
        raise SystemExit(0)
    main(sys.argv[1], gram="--no-gram" not in sys.argv, save=save)

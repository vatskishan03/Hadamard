"""Does an orthogonal design of order 4 or 8 admit a Goethals-Seidel-type circulant realization?

Entry (i,j) = sign * X_k^{s} R^{r}  (X^1 = transpose, R = back identity). For circulants X, Y:
    (X R^a)(Y R^b)^T = X Y^T            if a + b even   ("type 0")
                     = X Y R            if a + b odd    ("type 1")
Cross terms between rows i, i' cancel pairwise (the OD already supplies opposite signs) iff in the two
columns where the variable pair {k, l} meets, the two products are the same type and the same monomial.
All conditions are linear over GF(2) in the bits r_ij, s_ij."""
import itertools
import numpy as np

def od_matrices():
    Q = [["a","b","c","d"],["-b","a","-d","c"],["-c","d","a","-b"],["-d","-c","b","a"]]   # quaternion/Williamson
    O = ["a b c d e f g h", "-b a d -c f -e -h g", "-c -d a b g h -e -f", "-d c -b a h -g f -e",
         "-e -f -g -h a b c d", "-f e -h g -b a -d c", "-g h e -f -c d a -b", "-h -g f e -d -c b a"]
    O = [r.split() for r in O]
    return {"4": Q, "8": O}

def parse(M):
    n = len(M); var = {}; sign = np.zeros((n, n), int); v = np.zeros((n, n), int)
    for i in range(n):
        for j in range(n):
            t = M[i][j]; sg = -1 if t.startswith("-") else 1; name = t.lstrip("-")
            var.setdefault(name, len(var)); v[i, j] = var[name]; sign[i, j] = sg
    return v, sign

def check_od(v, sign):
    n = len(v); rng = np.random.default_rng(0); x = rng.standard_normal(n)
    A = sign * x[v]
    return np.allclose(A @ A.T, (x @ x) * np.eye(n))

def solve_gf2(rows, rhs, nvar):
    A = np.array(rows, dtype=np.uint8); b = np.array(rhs, dtype=np.uint8)
    M = np.concatenate([A, b[:, None]], axis=1); r = 0
    for c in range(nvar):
        p = next((i for i in range(r, len(M)) if M[i, c]), None)
        if p is None: continue
        M[[r, p]] = M[[p, r]]
        for i in range(len(M)):
            if i != r and M[i, c]: M[i] ^= M[r]
        r += 1
    inconsistent = any(not M[i, :nvar].any() and M[i, nvar] for i in range(len(M)))
    return not inconsistent, r

for name, M in od_matrices().items():
    v, sign = parse(M); n = len(v)
    assert check_od(v, sign), f"{name} is not an orthogonal design"
    R = lambda i, j: i * n + j               # r_ij bit index
    S = lambda i, j: n * n + i * n + j       # s_ij bit index
    rows, rhs = [], []
    def eq(idx, b):
        row = [0] * (2 * n * n)
        for t in idx: row[t] ^= 1
        rows.append(row); rhs.append(b)
    for i, ip in itertools.combinations(range(n), 2):
        meet = {}
        for j in range(n):
            meet.setdefault(frozenset((v[i, j], v[ip, j])), []).append(j)
        for key, cols in meet.items():
            assert len(cols) == 2, "each variable pair must meet in exactly two columns"
            j, jp = cols
            if v[i, j] != v[ip, jp]: raise ValueError("unexpected pattern")
            eq([R(i, j), R(ip, j), R(i, jp), R(ip, jp)], 0)                  # same type
            eq([S(i, j), S(ip, jp), R(i, j), R(ip, j)], 1)                   # same star on variable k
            eq([S(ip, j), S(i, jp), R(i, j), R(ip, j)], 1)                   # same star on variable l
    ok, rank = solve_gf2(rows, rhs, 2 * n * n)
    print(f"OD order {name}: {len(rows)} GF(2) equations in {2*n*n} bits, rank {rank}: "
          f"{'SOLVABLE - a GS-type circulant array exists' if ok else 'INCONSISTENT - no GS-type array'}")

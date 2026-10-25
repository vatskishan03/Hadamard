import numpy as np, sys
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form
def lat_index(K, ncols):
    D = (K[1:] - K[0])[:, :ncols].astype(object)
    # integer row reduction (HNF-like) using python ints
    A = [list(map(int, r)) for r in D]
    rows = []
    col = 0; r0 = 0
    A = [r for r in A if any(r)]
    piv_rows = []
    for c in range(ncols):
        # gcd-reduce column c among remaining rows
        while True:
            nz = [i for i in range(r0, len(A)) if A[i][c] != 0]
            if len(nz) <= 1: break
            nz.sort(key=lambda i: abs(A[i][c]))
            p = nz[0]
            for i in nz[1:]:
                q = A[i][c] // A[p][c]
                A[i] = [a - q * b for a, b in zip(A[i], A[p])]
        nz = [i for i in range(r0, len(A)) if A[i][c] != 0]
        if nz:
            A[r0], A[nz[0]] = A[nz[0]], A[r0]; r0 += 1
    diag = [abs(A[i][i]) if i < r0 else 0 for i in range(min(r0, ncols))]
    return r0, int(np.prod([float(x) for x in diag])), diag
if __name__ == "__main__":
    f, C, nc = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    K = np.fromfile(f, dtype=np.int16).reshape(-1, C)[:400].astype(int)
    print(lat_index(K, nc))

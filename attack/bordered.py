"""Bordered Goethals-Seidel arrays with Q = Z_3 and T = all characters of Q (border of size 12).

G is an abelian group of order n - 3 with a subgroup K of index 3; here G = Z_3 x K and K = 0 x K. Four +-1
sequences x_1..x_4 on G with

    sum_i PAF_i(d) = -12 for d in K \\ 0,     sum_i PAF_i(d) = 0 for d not in K,

together with a border (X, Rh, Ch) give a Hadamard matrix of order 4n:

    H = [ X  R ]      R[x, (i, g)] = Rh[x, (i, g + K)],   C[(i, g), x] = Ch[(i, g + K), x],
        [ C  M ]      M = the Goethals-Seidel array of x_1..x_4 over G (verify.build_bordered_gs's layout).

The border only depends on the compressions y_i(a) = sum_{k in K} x_i(a, k) (a in Z_3): with Mc the
Goethals-Seidel array of y_1..y_4 over Z_3, it needs X, Rh, Ch Hadamard of order 12 with X Ch^T + Rh Mc^T = 0
(for |K| > 18). The compressions satisfy sum_i PAF_{y_i} = 12 delta on Z_3; when |K| is even they are three
single +-2 entries in three different blocks and one zero block (4 x 27 x 8 = 864 patterns).
attack/border3.c solves the border for a given pattern (X fixed to Paley's matrix of order 12: all Hadamard matrices
of order 12 are equivalent, and (X, Rh, Ch) -> (P X Q, P Rh, Ch Q) preserves the conditions).

  python3 attack/bordered.py patterns            solve the border for all 864 compression patterns
  python3 attack/bordered.py validate N K1,K2..   exhaustive search on G = Z_3 x Z_K1 x ... (n = 3 + |G|) for the
                                                  first solutions with rows (0, 2, 2, 2); each is bordered,
                                                  assembled and checked (H H^T = 4N I)
"""
import itertools, subprocess, sys
from pathlib import Path as _P
import numpy as np
sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
import spence_lib as L


def gs_array(G, seqs):
    """Goethals-Seidel array over the abelian group G (layout of verify.build_bordered_gs)."""
    co = G.coords
    one = G.index(co[None, :, :] - co[:, None, :])
    two = G.index(co[:, None, :] + co[None, :, :])
    neg = G.index(-co)
    a, b, c, d = seqs
    A = a[one]; B, C, D = (x[two] for x in (b, c, d)); Br, Cr, Dr = (x[neg][two] for x in (b, c, d))
    return np.block([[A, B, C, D], [-B, A, Dr, -Cr], [-C, -Dr, A, Br], [-D, Cr, -Br, A]])


def paley12():
    p = 11; chi = [0] + [1 if pow(x, (p - 1) // 2, p) == 1 else -1 for x in range(1, p)]
    Qm = np.array([[chi[(j - i) % p] for j in range(p)] for i in range(p)])
    S = np.zeros((12, 12), dtype=np.int64); S[0, 1:] = 1; S[1:, 0] = -1; S[1:, 1:] = Qm
    return S + np.eye(12, dtype=np.int64)


def border(yc, nodes=2_000_000_000):
    """Border (X, Rh, Ch) for compressions yc (4 sequences on Z_3, entries 0 or +-2), or None."""
    Mc = gs_array(L.AbGroup([3]), [np.asarray(y, dtype=np.int64) for y in yc])
    X = paley12(); Wm = Mc // 2
    inp = '12\n' + ' '.join(map(str, X.ravel())) + '\n' + ' '.join(map(str, Wm.ravel())) + '\n'
    out = subprocess.run([W.bin('border3'), str(nodes)], input=inp, capture_output=True, text=True).stdout.split()
    if not out or out[0] != 'FOUND':
        return None
    Ch = np.array([[-1 if (int(r) >> i) & 1 else 1 for i in range(12)] for r in out[1:]], dtype=np.int64)
    R = -(X @ Ch.T @ Wm)
    assert np.all(np.abs(R) == 6)
    Rh = R // 6
    assert np.array_equal(Rh @ Rh.T, 12 * np.eye(12)) and np.array_equal(Ch @ Ch.T, 12 * np.eye(12))
    assert np.array_equal(X @ Ch.T + Rh @ Mc.T, np.zeros((12, 12)))
    return X, Rh, Ch


def assemble(G, seqs, X, Rh, Ch):
    """The bordered matrix of order 4|G| + 12 (G = Z_3 x K, coset index = first coordinate)."""
    v = G.v; cos = G.coords[:, 0] % 3
    M = gs_array(G, seqs)
    R = np.zeros((12, 4 * v), dtype=np.int64); C = np.zeros((4 * v, 12), dtype=np.int64)
    for i in range(4):
        R[:, i * v:(i + 1) * v] = Rh[:, 3 * i + cos]
        C[i * v:(i + 1) * v, :] = Ch[3 * i + cos, :]
    return np.block([[X, R], [C, M]])


def all_patterns():
    ok = 0; total = 0
    for zero in range(4):
        for pos in itertools.product(range(3), repeat=3):
            for sg in itertools.product((1, -1), repeat=3):
                yc = []; it = iter(zip(pos, sg))
                for i in range(4):
                    y = [0, 0, 0]
                    if i != zero:
                        p, s = next(it); y[p] = 2 * s
                    yc.append(y)
                total += 1; ok += border(yc) is not None
    print(f"{ok} of {total} compression patterns have a border")


def pool_vectors(st, r):
    O = st.O; out = []; sz = st.sizes.astype(np.int64)
    for start in range(0, 1 << O, 1 << 20):
        idx = np.arange(start, min(start + (1 << 20), 1 << O), dtype=np.int64)
        Y = (1 - 2 * ((idx[:, None] >> np.arange(O)) & 1)).astype(np.int8)
        out.append(Y[(Y.astype(np.int64) @ sz) == r])
    return np.concatenate(out)


def validate(n, kmods, maxsol=3):
    mods = [3] + kmods; G = L.AbGroup(mods); co = G.coords; v = G.v
    assert v == n - 3, "need n = 3 + |G|"
    gam = np.where(co[:, 0] == 0, 12, 0)
    st = L.Structure(G, [G.perm_from_matrix(np.eye(len(mods), dtype=np.int64))])
    tgt = np.array([-gam[st.reps[c]] for c in range(st.C)], dtype=np.int64)
    pools = {r: pool_vectors(st, r) for r in (0, 2)}
    vals = {r: np.rint(st.class_values(pools[r].astype(float))).astype(np.int64) for r in pools}
    A, B = vals[0], vals[2]
    pairs = {}
    for i in range(len(A)):
        for j in range(len(B)):
            pairs.setdefault((A[i] + B[j]).tobytes(), []).append((i, j))
    sols = []
    for k in range(len(B)):
        for l in range(len(B)):
            for (i, j) in pairs.get((tgt - B[k] - B[l]).tobytes(), []):
                sols.append((i, j, k, l))
                if len(sols) >= maxsol: break
            if len(sols) >= maxsol: break
        if len(sols) >= maxsol: break
    for (i, j, k, l) in sols:
        seqs = [pools[0][i], pools[2][j], pools[2][k], pools[2][l]]
        seqs = [s.astype(np.int64)[st.orb] for s in seqs]
        for d in range(1, v):
            sh = G.index(co + co[d])
            assert sum(int(x @ x[sh]) for x in seqs) == -gam[d]
        yc = [[int(s[co[:, 0] == a].sum()) for a in range(3)] for s in seqs]
        b = border(yc)
        if b is None:
            print(f"  compressions {yc}: no border"); continue
        H = assemble(G, seqs, *b)
        good = np.array_equal(H @ H.T, 4 * n * np.eye(4 * n, dtype=np.int64))
        print(f"n = {n}, G = Z_3 x {' x '.join(f'Z_{m}' for m in kmods)}: compressions {yc}, "
              f"H of order {4 * n} with H H^T = {4 * n} I: {good}")


if __name__ == '__main__':
    if sys.argv[1] == 'patterns':
        all_patterns()
    elif sys.argv[1] == 'validate':
        validate(int(sys.argv[2]), [int(x) for x in sys.argv[3].split(',')])

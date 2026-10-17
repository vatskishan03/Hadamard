"""Skew-Hadamard matrices of order h m from Wallis's product of a skew-Hadamard matrix of order h = 4v and an amicable
pair of order m (J. Wallis, J. Combin. Theory Ser. A 21 (1976); Seberry-Yamada 1992, SXI): K = I (x) A + S (x) M.

    python3 search/skew_wallis.py paley3 v h      skew Goethals-Seidel quadruples of prime length v = 3 (mod 4):
                                                  the Paley sequence first, three blocks invariant under the subgroup
                                                  of order h of the squares (exhaustive over unions of orbits)
    python3 search/skew_wallis.py certify order   writes skew/<order>/certificate.json and sequences.txt
                                                  (orders 9424, 9776, 13072)

The amicable pairs of orders 52 and 76 are those of skew/2652 and skew/5700 (skew-propus quadruples of lengths 13
and 19). verify/verify.py (build_skew_amicable_product) rebuilds and checks the matrices.
"""
import hashlib
import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "verify"))
import verify as V  # noqa: E402


def squares(v):
    return sorted({x * x % v for x in range(1, v)})


def subgroup(v, h):
    """The subgroup of order h of Z_v^* (v prime, h | v - 1)."""
    g = next(x for x in range(2, v) if len({pow(x, k, v) for k in range(v - 1)}) == v - 1)
    return sorted({pow(g, (v - 1) // h * k, v) for k in range(h)})


def paley3(v, h, limit=2):
    """Three +-1 sequences of length v, unions of orbits of the subgroup H of order h of the squares, with
    PAF_1 + PAF_2 + PAF_3 = 1 off zero: with the Paley sequence (PAF -1 off zero) a skew Goethals-Seidel quadruple
    (the difference families of Djokovic-Kotsireas, arXiv:1511.08734). Returns the first `limit` triples of sets."""
    H = subgroup(v, h)
    if not set(H) <= set(squares(v)):
        raise ValueError("H must lie in the squares")
    orbits, seen = [[0]], {0}
    for x in range(1, v):
        if x not in seen:
            o = sorted({x * y % v for y in H})
            seen |= set(o)
            orbits.append(o)
    reps = [o[0] for o in orbits[1:]]                     # the PAF is H-invariant: check it on orbit representatives
    nb = len(orbits)
    X = np.ones((2 ** nb, v), dtype=np.int64)
    for m in range(2 ** nb):
        for b in range(nb):
            if m >> b & 1:
                X[m, orbits[b]] = -1
    P = np.rint(np.fft.ifft(np.abs(np.fft.fft(X, axis=1)) ** 2, axis=1).real).astype(np.int64)[:, reps]
    rs = np.abs(X.sum(1))
    by_row = {r: np.flatnonzero(rs == r) for r in sorted(set(rs.tolist()))}
    out = []
    for r1 in by_row:
        for r2 in by_row:
            r3sq = 4 * v - 1 - r1 * r1 - r2 * r2
            r3 = int(round(max(r3sq, 0) ** 0.5))
            if r2 < r1 or r3 < r2 or r3 * r3 != r3sq or r3 not in by_row:
                continue
            table = {}
            for k in by_row[r3]:
                table.setdefault(P[k].tobytes(), k)
            L2 = by_row[r2]
            for i in by_row[r1]:
                for j, row in zip(L2, 1 - P[i][None, :] - P[L2]):
                    k = table.get(row.tobytes())
                    if k is not None:
                        out.append([sorted(np.flatnonzero(X[t] < 0).tolist()) for t in (i, j, k)])
                        if len(out) == limit:
                            return H, out
    return H, out


def orbit_reps(s, H, v):
    s, reps = set(s), []
    while s:
        r = min(s)
        o = {r * y % v for y in H}
        if not o <= s:
            raise ValueError("set is not H-invariant")
        reps.append(r)
        s -= o
    return reps


def factor_cert(v, h, sets):
    H = subgroup(v, h) if h > 1 else [1]
    return {"v": v, "multiplier": H[1] if h > 1 else 1, "subgroup_order": h,
            "blocks": [orbit_reps(squares(v), H, v)] + [orbit_reps(s, H, v) for s in sets],
            "rowsums": [v - 2 * len(squares(v))] + [v - 2 * len(s) for s in sets]}


# Djokovic, arXiv:0704.0640 (Int. Math. Forum 3 (2008)), Proposition 2.1: X1, X2, X3 in Z_47, with X0 the squares
DJOKOVIC_47 = [[2, 3, 5, 6, 7, 9, 10, 11, 12, 13, 14, 17, 18, 19, 20, 21, 22, 25, 27, 30, 31, 33, 35, 37, 38, 39, 40,
                42, 43, 44],
               [1, 3, 6, 7, 8, 11, 13, 14, 15, 19, 20, 21, 24, 27, 30, 33, 39, 41, 43, 44, 45, 46],
               [3, 6, 8, 10, 11, 12, 14, 20, 21, 23, 24, 25, 26, 27, 30, 31, 32, 34, 35, 41, 42, 45]]

ORDERS = {9424: (31, 5, "5700"), 13072: (43, 7, "5700"), 9776: (47, 1, "2652")}


def certify(order):
    v, h, amicable_dir = ORDERS[order]
    sets = DJOKOVIC_47 if v == 47 else paley3(v, h, 1)[1][0]
    am = json.loads((ROOT / "skew" / amicable_dir / "certificate.json").read_text())
    cert = {"order": order, "construction": "skew-amicable-product", "skew": True, "factor": factor_cert(v, h, sets),
            "amicable": {"n": am["n"], "blocks": am["blocks"], "rowsums": am["rowsums"]}}
    K, seqs = V.build_skew_amicable_product(cert)
    if not (K.shape[0] == order and V.is_pm1(K) and V.gram_is_scalar(K, order)
            and np.array_equal(K + K.T, 2 * np.eye(order, dtype=K.dtype))):
        raise ValueError("not a skew-Hadamard matrix")
    sha = hashlib.sha256()
    for chunk in V.csv_chunks(K):
        sha.update(chunk)
    cert["sha256"] = sha.hexdigest()
    d = ROOT / "skew" / str(order)
    d.mkdir(exist_ok=True)
    (d / "certificate.json").write_text(json.dumps(cert) + "\n")
    (d / "sequences.txt").write_text("".join("".join("+" if t > 0 else "-" for t in x) + "\n" for x in seqs))
    print(order, "ok", cert["sha256"])


if __name__ == "__main__":
    if sys.argv[1] == "paley3":
        v, h = int(sys.argv[2]), int(sys.argv[3])
        H, sols = paley3(v, h)
        for s in sols:
            print("v", v, "H", H, "row sums", [v - 2 * len(x) for x in s], "sets", s)
        if not sols:
            print("none")
    elif sys.argv[1] == "certify":
        for o in sys.argv[2:]:
            certify(int(o))

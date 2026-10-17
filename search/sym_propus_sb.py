#!/usr/bin/env python3
"""Symmetric Hadamard matrices of order 4(q + 2) by Seberry and Balonin's Corollary 1 (Australas. J. Combin. 69 (2017);
the census's Seberry-Balonin series), for a prime q = 1 (mod 4) and the core Q of a symmetric conference matrix of
order (q + 1)/2 + 1: Seberry-Yamada's Lemma 26 with U = (I, 0, Q, Q), V = (M, N, I, I), M and N Yamada's symmetric
circulants with M^2 + N^2 = qI, then the propus array (verify.py construction `symmetric-propus-sb`).

    python3 search/sym_propus_sb.py certify v CONF.npy   writes symmetric/<4v>/certificate.json from a symmetric
                                                         conference matrix of order (v - 1)/2 + 1 stored by numpy
    python3 search/sym_propus_sb.py mathon t OUT.npy     Mathon's conference matrix of order (4t - 1)^2 (4t + 1) + 1
                                                         from SageMath (passagemath-graphs), first row and column 1
    python3 search/sym_propus_sb.py skewsq m OUT.npy     the conference matrix of order (m - 1)^2 + 1 from a
                                                         skew-Hadamard matrix of order m

v = 883 (q = 881) uses Mathon's matrix of order 442 (t = 2), v = 451 (q = 449) the matrix of order 226 (m = 16).
"""
import base64, hashlib, json, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path[:0] = [str(ROOT / "verify"), str(ROOT / "search")]
import verify  # noqa: E402
from make_yamada import primitive  # noqa: E402


def mathon(t, out):
    import sage.all__sagemath_graphs  # noqa: F401  (initialises SageMath)
    from sage.graphs.generators.families import MathonPseudocyclicStronglyRegularGraph
    A = np.array(MathonPseudocyclicStronglyRegularGraph(t).seidel_adjacency_matrix(), dtype=np.int64)
    n = A.shape[0]
    C = np.zeros((n + 1, n + 1), dtype=np.int64)
    C[0, 1:] = C[1:, 0] = 1
    C[1:, 1:] = A
    assert np.array_equal(C, C.T) and np.array_equal(C @ C, n * np.eye(n + 1, dtype=np.int64))
    np.save(out, C)
    print(f"symmetric conference matrix of order {n + 1} saved to {out}")


def skewsq(m, out):
    """the symmetric conference matrix of order (m - 1)^2 + 1 from a skew-Hadamard matrix of order m (Brouwer's
    "skewhad^2"): with the skew core S of order k = m - 1 (S^T = -S, S J = 0, S S^T = kI - J), the core is
    Q = S (x) S + I (x) J - J (x) I. m = 2^a * (p + 1) with p = 3 (mod 4) a prime: Paley's skew matrix, doubled."""
    import math
    p, a = m, 0
    while p % 2 == 0 and not (all((p - 1) % d for d in range(2, math.isqrt(p - 1) + 1)) and (p - 1) % 4 == 3):
        p, a = p // 2, a + 1
    p -= 1
    chi = np.array([0] + [1 if pow(x, (p - 1) // 2, p) == 1 else -1 for x in range(1, p)], dtype=np.int64)
    H = np.ones((p + 1, p + 1), dtype=np.int64)
    H[1:, 0] = -1
    H[1:, 1:] = np.eye(p, dtype=np.int64) + chi[(np.arange(p)[None, :] - np.arange(p)[:, None]) % p]
    for _ in range(a):                                       # [[H, H], [-H^T, H^T]] is skew-Hadamard again
        H = np.block([[H, H], [-H.T, H.T]])
    k = m - 1
    assert H.shape == (m, m) and np.array_equal(H + H.T, 2 * np.eye(m, dtype=np.int64))
    s = H[0].copy()                                          # normalise: first row all 1, keeping H + H^T = 2I
    H = H * s[None, :] * s[:, None]
    S = H[1:, 1:] - np.eye(k, dtype=np.int64)
    I, J = np.eye(k, dtype=np.int64), np.ones((k, k), dtype=np.int64)
    Q = np.kron(S, S) + np.kron(I, J) - np.kron(J, I)
    n = k * k
    C = np.zeros((n + 1, n + 1), dtype=np.int64)
    C[0, 1:] = C[1:, 0] = 1
    C[1:, 1:] = Q
    assert np.array_equal(C, C.T) and np.array_equal(C @ C, n * np.eye(n + 1, dtype=np.int64))
    np.save(out, C)
    print(f"symmetric conference matrix of order {n + 1} saved to {out}")


def certify(v, conf):
    q = v - 2
    n = (q + 1) // 2
    C = np.load(conf)
    if C.shape != (n + 1, n + 1):
        raise SystemExit(f"need a conference matrix of order {n + 1}")
    s = np.where(C[0, 1:] < 0, -1, 1)              # normalise the first row and column to 1
    C = C * np.concatenate(([1], s))[:, None] * np.concatenate(([1], s))[None, :]
    Q = C[1:, 1:]
    bits = (Q[np.triu_indices(n, 1)] > 0).astype(np.uint8)
    core = {"order": n, "upper": base64.b64encode(np.packbits(bits).tobytes()).decode(),
            "sha256": hashlib.sha256(Q.astype(np.int8).tobytes()).hexdigest()}
    r = next(x for x in range(2, q) if pow(x, (q - 1) // 2, q) == q - 1)
    cert = {"order": 4 * v, "construction": "symmetric-propus-sb", "symmetric": True, "q": q, "r": r,
            "xi": list(primitive(q, r)), "core": core}
    H, _ = verify.BUILDERS[cert["construction"]](cert)
    m = H.shape[0]
    assert m == 4 * v and verify.is_pm1(H) and verify.gram_is_scalar(H, m) and (H == H.T).all(), "not symmetric Hadamard"
    sha = hashlib.sha256()
    for chunk in verify.csv_chunks(H):
        sha.update(chunk)
    cert["sha256"] = sha.hexdigest()
    d = ROOT / "symmetric" / str(4 * v)
    d.mkdir(exist_ok=True)
    (d / "certificate.json").write_text(json.dumps(cert, indent=1) + "\n")
    print(f"symmetric {4 * v} = 4 * {v}: H = H^T, H H^T = {m} I verified; sha256 {cert['sha256']}", flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "mathon":
        mathon(int(sys.argv[2]), sys.argv[3])
    elif sys.argv[1] == "skewsq":
        skewsq(int(sys.argv[2]), sys.argv[3])
    elif sys.argv[1] == "certify":
        certify(int(sys.argv[2]), sys.argv[3])

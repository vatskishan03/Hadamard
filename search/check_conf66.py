#!/usr/bin/env python3
"""Check the two symmetric conference matrices of order 66 in search/data/conference66_*.txt and the 16 strongly
regular graphs srg(65, 32, 15, 16) in search/data/srg65_32_15_16_new16.g6 (search/NOTES.md §62).

    python3 search/check_conf66.py

For each matrix C: C = C^T, zero diagonal, +-1 off it, C^2 = 65 I. Its descendants (switch C or -C so that row b is
all +1 off the diagonal, delete b, A = (J - I - S)/2) are srg(65, 32, 15, 16); the script counts the distinct ones
(up to isomorphism if pynauty is installed, else up to equality) and checks that the g6 file lists exactly them.
"""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "search" / "data"


def g6_decode(s):
    b = [ord(c) - 63 for c in s.strip()]
    if b[0] == 63:
        n, rest = (b[1] << 12) | (b[2] << 6) | b[3], b[4:]
    else:
        n, rest = b[0], b[1:]
    bits = [(x >> k) & 1 for x in rest for k in range(5, -1, -1)]
    A = np.zeros((n, n), dtype=np.int64)
    t = 0
    for j in range(1, n):
        for i in range(j):
            A[i, j] = A[j, i] = bits[t]
            t += 1
    return A


def is_srg(A):
    n = len(A)
    return (A == A.T).all() and not np.diag(A).any() and (A.sum(1) == 32).all() and \
        (A @ A + A == 16 * np.eye(n, dtype=np.int64) + 16).all()


def key(A):
    try:
        import pynauty
        return pynauty.certificate(pynauty.Graph(len(A), adjacency_dict={i: list(np.flatnonzero(A[i])) for i in range(len(A))}))
    except ImportError:
        return A.tobytes()


def main():
    graphs = [g6_decode(l) for l in (DATA / "srg65_32_15_16_new16.g6").read_text().splitlines() if l.strip()]
    assert len(graphs) == 16 and all(is_srg(A) for A in graphs)
    listed = {key(A) for A in graphs}
    found = set()
    for f in sorted(DATA.glob("conference66_*.txt")):
        C = np.array([[{"0": 0, "+": 1, "-": -1}[c] for c in l] for l in f.read_text().split()], dtype=np.int64)
        n = len(C)
        assert n == 66 and (C == C.T).all() and not np.diag(C).any() and (np.abs(C) + np.eye(n, dtype=np.int64) == 1).all()
        assert (C @ C == 65 * np.eye(n, dtype=np.int64)).all()
        desc = set()
        for s in (1, -1):
            for b in range(n):
                d = s * C[b].copy(); d[b] = 1
                S = (d[:, None] * (s * C) * d[None, :])
                idx = [i for i in range(n) if i != b]
                S = S[np.ix_(idx, idx)]
                A = (1 - S) // 2
                np.fill_diagonal(A, 0)
                assert is_srg(A)
                desc.add(key(A))
        print(f"{f.name}: symmetric conference matrix of order 66, {len(desc)} distinct descendants srg(65,32,15,16)")
        found |= desc
    print("descendants listed in srg65_32_15_16_new16.g6:", found == listed)


if __name__ == "__main__":
    main()

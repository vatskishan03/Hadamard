#!/usr/bin/env python3
"""Exhaustive search for pairwise amicable T-matrices of small odd order t (NOTES.md §106).

T-matrices T_1, ..., T_4 of order t are circulant over {0, 1, -1} with disjoint supports covering Z_t and
sum T_k T_k^T = tI. Cooper-Wallis multiplication of an arbitrary Goethals-Seidel quadruple (not Williamson-type) would
need T_j T_k^T = T_k T_j^T for all j != k. Every assignment of a matrix and a sign to each position is tried (position
0 in T_1 with sign +, matrices labelled in order of first use). Output: the number of solutions; none for t = 3, 5, 7.

    python3 search/amicable_tmatrices.py 3 5 7
"""
import itertools
import sys

import numpy as np


def circ(v):
    t = len(v)
    return np.array([[v[(j - i) % t] for j in range(t)] for i in range(t)])


def count(t):
    found = 0
    for lab in itertools.product(range(4), repeat=t):
        first = []
        for x in lab:
            if x not in first:
                first.append(x)
        if lab[0] != 0 or first != sorted(first):
            continue
        for sg in itertools.product((1, -1), repeat=t - 1):
            T = [np.zeros(t, dtype=np.int64) for _ in range(4)]
            for x, s in enumerate((1,) + sg):
                T[lab[x]][x] = s
            C = [circ(v) for v in T]
            if not np.array_equal(sum(c @ c.T for c in C), t * np.eye(t, dtype=np.int64)):
                continue
            if all(np.array_equal(C[j] @ C[k].T, C[k] @ C[j].T) for j in range(4) for k in range(j + 1, 4)):
                found += 1
    return found


if __name__ == "__main__":
    for t in map(int, sys.argv[1:] or ["3", "5", "7"]):
        print(f"t = {t}: {count(t)} sets of pairwise amicable T-matrices")

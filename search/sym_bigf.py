"""Symmetric Hadamard matrices of order 4q with an involution of 2q - 2 fixed points (NOTES.md §54).

A symmetric H of order 4q with a permutation automorphism σ of order 2, f = 2q - 2 fixed points and t = q + 1
two-cycles is, in the basis (fixed points, first points of the pairs, second points),

    H = [[F, B, B], [Bᵀ, P + V, P - V], [Bᵀ, P - V, P + V]],

with V a symmetric weighing matrix of order q + 1 and weight q (V² = qI), P a symmetric signed permutation matrix on
the zero pattern of V (an involution π), and the Q-part

    BᵀB = f I,   F B = -2 B P,   F² + 2 B Bᵀ = 4q I,   F symmetric ±1.

For q ≡ 3 (mod 4) a prime power, Paley's skew core K and R: x ↦ -x (with -1 at ∞) give V = KR, with π fixing 0 and
∞; only the Q-part would then depend on q. For q = 7, B is 12 × 8 with orthogonal columns, so it is 8 columns of the
unique H(12) (every such B extends; M₁₂ is 5-transitive on the columns), and each row of F lies in the affine set
{x ∈ {±1}¹²: Bᵀx = -2(BP)ᵢ}. This program enumerates those sets for every involution π on the 8 columns and every
sign pattern of P, and backtracks on symmetry and the Gram condition.

    python3 search/sym_bigf.py 3      # control: Q-parts exist (order 12)
    python3 search/sym_bigf.py 7      # all π: Q-parts only for π with 4 fixed points, and no V for those
"""
import itertools
import sys

import numpy as np


def paley1(p):
    """Skew Hadamard matrix of order p + 1 (p ≡ 3 mod 4 prime)."""
    S = np.eye(p + 1, dtype=int)
    S[0, 1:] = 1
    S[1:, 0] = -1
    for x in range(p):
        for y in range(p):
            if x != y:
                S[1 + x, 1 + y] = 1 if pow((y - x) % p, (p - 1) // 2, p) == 1 else -1
    return S


def row_candidates(B, target):
    """All x ∈ {±1}^f with Bᵀx = target (meet in the middle over the rows)."""
    f = B.shape[0]
    h = f // 2
    left = {}
    for bits in itertools.product((1, -1), repeat=h):
        left.setdefault(tuple(np.array(bits) @ B[:h]), []).append(bits)
    out = []
    for bits in itertools.product((1, -1), repeat=f - h):
        for lb in left.get(tuple(target - np.array(bits) @ B[h:]), ()):
            out.append(np.array(lb + bits))
    return out


def q_parts(B, P, q):
    """All symmetric ±1 F with F B = -2 B P and F² + 2 B Bᵀ = 4q I."""
    f = B.shape[0]
    G = 4 * q * np.eye(f, dtype=int) - 2 * B @ B.T
    T = -2 * B @ P
    cands = [row_candidates(B, T[i]) for i in range(f)]
    if min(len(c) for c in cands) == 0:
        return []
    order = sorted(range(f), key=lambda i: len(cands[i]))
    F = np.zeros((f, f), dtype=int)
    done, sols = [], []

    def bt(k):
        if k == f:
            sols.append(F.copy())
            return
        i = order[k]
        for x in cands[i]:
            if all(x[j] == F[j, i] and x @ F[j] == G[i, j] for j in done):
                F[i] = x
                done.append(i)
                bt(k + 1)
                done.pop()
        F[i] = 0

    bt(0)
    return sols


def matchings(L):
    if not L:
        yield []
        return
    for i in range(1, len(L)):
        for m in matchings(L[1:i] + L[i + 1:]):
            yield [(L[0], L[i])] + m


def weighing_on_pattern(zero, q):
    """A symmetric ±1 matrix V, 0 on `zero`, with V² = qI (CP-SAT), or None."""
    from ortools.sat.python import cp_model
    n = len(zero)
    m = cp_model.CpModel()
    x = {(i, j): m.NewBoolVar('') for i in range(n) for j in range(i, n) if not zero[i, j]}

    def v(i, j):
        return x.get((min(i, j), max(i, j)))

    for i in range(n):
        for l in range(i + 1, n):
            terms = []
            for k in range(n):
                a, b = v(i, k), v(l, k)
                if a is not None and b is not None:
                    p = m.NewBoolVar('')
                    m.Add(a == b).OnlyEnforceIf(p)
                    m.Add(a != b).OnlyEnforceIf(p.Not())
                    terms.append(p)
            m.Add(2 * sum(terms) == len(terms))
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = 600
    r = s.Solve(m)
    if r == cp_model.INFEASIBLE:
        return None
    if r not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        raise RuntimeError('undecided')
    V = np.zeros((n, n), dtype=int)
    for (i, j), b in x.items():
        V[i, j] = V[j, i] = 1 if s.Value(b) else -1
    return V


def main(q):
    f, n = 2 * q - 2, q + 1
    B = paley1(f - 1)[:, :n]
    for k in range(n % 2, n + 1, 2):
        tried = found = complete = 0
        for fixed in itertools.combinations(range(n), k):
            rest = [c for c in range(n) if c not in fixed]
            for signs in itertools.product((1, -1), repeat=k):
                if k and signs[0] == -1:
                    continue  # P ↦ -P, F ↦ -F
                for mt in matchings(rest):
                    P = np.zeros((n, n), dtype=int)
                    for c, s in zip(fixed, signs):
                        P[c, c] = s
                    for a, b in mt:
                        P[a, b] = P[b, a] = 1  # pair signs are absorbed by column signs of B
                    tried += 1
                    sols = q_parts(B, P, q)
                    if not sols:
                        continue
                    found += 1
                    V = weighing_on_pattern(P != 0, q) if k < n else None
                    if V is None:
                        continue
                    F = sols[0]
                    H = np.block([[F, B, B], [B.T, P + V, P - V], [B.T, P - V, P + V]])
                    assert (H == H.T).all() and (abs(H) == 1).all() and (H @ H == 4 * q * np.eye(4 * q)).all()
                    complete += 1
        print(f'q = {q}, π with {k} fixed points: {tried} labellings, {found} with a Q-part, '
              f'{complete} completed to a symmetric H({4 * q})', flush=True)


if __name__ == '__main__':
    main(int(sys.argv[1]))

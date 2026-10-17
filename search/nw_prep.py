#!/usr/bin/env python3
"""Tables for nw_exh.c: near-Williamson quadruples on Z_v with a multiplier group H (|H| odd, -1 not in H).

A near-Williamson quadruple is (A, B, C, D), circulant +-1 of order v, B, C, D symmetric, PAF_A + PAF_B + PAF_C +
PAF_D = 4v delta. With all four constant on the orbits of H (B, C, D then on the orbits of L = <H, -1>) the condition
holds on the L-orbits of Z_v \\ {0} ("classes").

nw_exh.c splits A on each L-orbit O u -O into a symmetric part (A equal on O and -O) and an antisymmetric part
(A opposite on O and -O). With A_s, A_a these two ternary sequences, PAF_A = PAF_{A_s} + PAF_{A_a} exactly (the
cross-correlation of a symmetric and an antisymmetric sequence is odd, so it cancels in the sum), and likewise
|A^(t)|^2 = A_s^(t)^2 + |A_a^(t)|^2.

Output (little-endian int32 unless noted), read by nw_exh.c:
    v, C, O (= number of L-orbits including {0}), rows rA rB rC rD
    sizesL[O]
    QL[C][O][O]        PAF of a symmetric L-invariant sequence y: sum_{i,j} y_i y_j QL[c][i][j]
    QA[C][O-1][O-1]    PAF of the antisymmetric part x (x_k on O_k+, -x_k on O_k-): sum_{k,l} x_k x_l QA[c][k][l]
    GL[C][O]  (float64) symmetric Fourier: y^(t_c) = sum_i y_i GL[c][i]
    GA[C][O-1] (float64) antisymmetric Fourier: |A_a^(t_c)| = |sum_k x_k GA[c][k]|  (GA = sum over O_k+ minus O_k- of sin)
    orbL[v], plusorb[v] (int32): L-orbit of each g, and +1/-1 for g in O_k+ / O_k- (0 for g = 0)

    python3 search/nw_prep.py v g [g2 ...] --rows rA rB rC rD --out FILE
"""
import argparse, struct, sys
import numpy as np


def closure(gens, m):
    S, todo = {1}, [1]
    while todo:
        x = todo.pop()
        for g in gens:
            y = x * g % m
            if y not in S:
                S.add(y); todo.append(y)
    return sorted(S)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("v", type=int); ap.add_argument("gens", type=int, nargs="+")
    ap.add_argument("--rows", type=int, nargs=4, required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    v = a.v
    H = closure(a.gens, v)
    if (v - 1) in H:
        sys.exit("-1 in H")
    orbL = -np.ones(v, np.int64); orbL[0] = 0
    plus = np.zeros(v, np.int64)
    reps = []
    for g in range(1, v):
        if orbL[g] < 0:
            k = len(reps) + 1
            reps.append(g)
            for h in H:
                orbL[g * h % v] = k; plus[g * h % v] = 1
                orbL[(-g * h) % v] = k; plus[(-g * h) % v] = -1
    O = len(reps) + 1
    C = O - 1
    sizesL = np.bincount(orbL, minlength=O)
    t = np.arange(v)
    QL = np.zeros((C, O, O), np.int64)
    QA = np.zeros((C, O - 1, O - 1), np.int64)
    nz = t != 0
    for c, s in enumerate(reps):
        u = (t + s) % v
        np.add.at(QL[c], (orbL, orbL[u]), 1)
        m = nz & (u != 0)
        np.add.at(QA[c], (orbL[m] - 1, orbL[u[m]] - 1), plus[m] * plus[u[m]])
    # checks against direct PAFs
    rng = np.random.default_rng(1)
    for _ in range(5):
        y = rng.choice([-1, 1], O)
        x = rng.choice([-1, 1], O - 1)
        Ys = y[orbL]
        Xa = np.zeros(v, np.int64); Xa[nz] = x[orbL[nz] - 1] * plus[nz]
        for c, s in enumerate(reps):
            assert int(Ys @ np.roll(Ys, -s)) == int(y @ QL[c] @ y)
            assert int(Xa @ np.roll(Xa, -s)) == int(x @ QA[c] @ x)
    ang = 2 * np.pi * np.outer(reps, t) / v           # characters t_c in the same class representatives
    GL = np.zeros((C, O)); GA = np.zeros((C, O - 1))
    for c in range(C):
        np.add.at(GL[c], orbL, np.cos(ang[c]))
        np.add.at(GA[c], orbL[nz] - 1, plus[nz] * np.sin(ang[c][nz]))
    # Fourier checks: |sum y_g zeta^{t g}| and the antisymmetric part
    for _ in range(3):
        y = rng.choice([-1, 1], O); x = rng.choice([-1, 1], O - 1)
        Ys = y[orbL]; Xa = np.zeros(v); Xa[nz] = x[orbL[nz] - 1] * plus[nz]
        F = np.fft.fft(Ys); Fa = np.fft.fft(Xa)
        for c, s in enumerate(reps):
            # numpy fft uses exp(-2 pi i t g / v); magnitudes agree
            assert abs(abs(F[s]) - abs(y @ GL[c])) < 1e-6
            assert abs(abs(Fa[s]) - abs(x @ GA[c])) < 1e-6
    rA, rB, rC, rD = a.rows
    if rA * rA + rB * rB + rC * rC + rD * rD != 4 * v:
        sys.exit("row sums do not satisfy a^2 + b^2 + c^2 + d^2 = 4v")
    with open(a.out, "wb") as f:
        f.write(struct.pack("<8i", v, C, O, rA, rB, rC, rD, len(H)))
        f.write(sizesL.astype("<i4").tobytes())
        f.write(QL.astype("<i4").tobytes())
        f.write(QA.astype("<i4").tobytes())
        f.write(GL.astype("<f8").tobytes())
        f.write(GA.astype("<f8").tobytes())
        f.write(orbL.astype("<i4").tobytes())
        f.write(plus.astype("<i4").tobytes())
    print(f"v={v} |H|={len(H)} C={C} O={O} rows={a.rows} -> {a.out}")


if __name__ == "__main__":
    main()

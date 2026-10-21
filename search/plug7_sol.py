#!/usr/bin/env python3
"""A SOL line of plug7_l2 (four blocks, each p columns 'X0:E') -> the four sequences of the zero-diagonal Szollosi array
over Z_{7p} (index g with g mod 7 the F_7 coordinate, g mod p the Z_p coordinate), checked: sum PAF = (28p - 2) delta,
A, B symmetric with zeros exactly at 0, all four invariant under H = squares of F_7.
X0 encodes x0 and the defect sign e (search/NOTES.md section 82): +-1 no defect; +-7: x0 = e = +-1; +-5: x0 = -e = -+1;
+-6: x0 = 0, e = +-1. With t = E, the block is x0 on {0}, e + t on the squares (1, 2, 4), e - t on the non-squares.
usage: plug7_sol.py SOLFILE p [line] [out.npy]   (prints the four sequences as +/-/0 strings)"""
import sys
import numpy as np

SQ = (1, 2, 4)
DEC = {1: (1, 0), -1: (-1, 0), 7: (1, 1), -7: (-1, -1), 5: (-1, 1), -5: (1, -1), 6: (0, 1), -6: (0, -1)}


def blocks(line, p):
    out = []
    for part in line.split("|")[1:]:
        cols = [tuple(map(int, c.split(":"))) for c in part.split()]
        assert len(cols) == p
        x0 = np.array([DEC[a][0] for a, _ in cols]); e = np.array([DEC[a][1] for a, _ in cols]); t = np.array([b for _, b in cols])
        x = np.zeros(7 * p, dtype=np.int64)
        for g in range(7 * p):
            y, k = g % 7, g % p
            x[g] = x0[k] if y == 0 else (e[k] + t[k] if y in SQ else e[k] - t[k])
        out.append(x)
    return np.array(out)


def main():
    fn, p = sys.argv[1], int(sys.argv[2]); which = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    line = [l for l in open(fn) if l.startswith("SOL")][which]
    S = blocks(line, p); m = 7 * p
    assert set(np.unique(S[2:])) <= {-1, 1} and all(S[k][0] == 0 and (S[k][1:] != 0).all() for k in (0, 1))
    assert all(np.array_equal(S[k][(-np.arange(m)) % m], S[k]) for k in (0, 1))
    g = next(x for x in range(2, m) if x % p == 1 and x % 7 == 2)
    assert all(np.array_equal(S[k][(g * np.arange(m)) % m], S[k]) for k in range(4))
    F = np.fft.fft(S.astype(float), axis=1)
    paf = np.rint(np.fft.ifft((np.abs(F) ** 2).sum(0)).real).astype(np.int64)
    assert paf[0] == 4 * m - 2 and not paf[1:].any()
    print(f"m = {m}: sum PAF = {4 * m - 2} delta, H = <{g}>, row sums {[int(v) for v in S.sum(1)]}")
    for x in S:
        print("".join({1: "+", -1: "-", 0: "0"}[int(v)] for v in x))
    if len(sys.argv) > 4:
        np.save(sys.argv[4], S)


if __name__ == "__main__":
    main()

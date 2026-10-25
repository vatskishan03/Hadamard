#!/usr/bin/env python3
"""Assemble the symmetric conference matrix of order 118 from an exact cover found by ./xc (see README.md).

    python3 build.py xc.out          reads insts.bin and the rows of the first SOLUTION in xc.out, builds C, checks
                                     C = C^T, zero diagonal, C^2 = 117 I exactly, normalises C (the extra point first,
                                     first row and column +) and writes conference118_<sha256 prefix>.txt

Points: 117 = 13 * 9 points (i, x) with i in Z_13 and x in AG(2, 3), plus one point. A row of insts.bin is a K_4
component: four points of Z_13 (its six pairs), with a sign D = +-1 at each (a "flag"). Each point of Z_13 lies in four
components, two with D = +1 (classes A1, A2) and two with D = -1 (B1, B2); each class is a direction class of AG(2, 3)
(A1 = (0, 1), A2 = (1, 0), B1 = (1, 1), B2 = (1, 2)). Diagonal blocks are J - I - 2R with R the "rook" graph of the two
B directions; the block of a pair (i, j) in a component is -chi(i - j) (J - 2M), M the hyperplane matrix
{(x, y): <a, x> - <b, y> = c} of the two flags' directions a, b, and c read off the cube-root phases that make the
component's 4 x 4 matrix W = D + 2 Phi satisfy W^2 = 13 I over Z[w]."""
import hashlib
import itertools
import struct
import sys

import numpy as np

P = 13
QR = {x * x % P for x in range(1, P)}
chi = lambda x: 0 if x % P == 0 else (1 if x % P in QR else -1)
PAIRS = [(i, j) for i in range(P) for j in range(i + 1, P)]
W3 = np.exp(2j * np.pi / 3)
PTS9 = [(a, b) for a in range(3) for b in range(3)]
J9, I9 = np.ones((9, 9), np.int64), np.eye(9, dtype=np.int64)
CLASSREP = {"A1": (0, 1), "A2": (1, 0), "B1": (1, 1), "B2": (1, 2)}


def read_rows(insts, xcout, types=1):
    data = open(insts, "rb").read()
    off, rows = 0, []
    for t in range(4):
        n = struct.unpack_from("<I", data, off)[0]
        off += 4
        for _ in range(n):
            lo, hi, us = struct.unpack_from("<QQQ", data, off)
            off += 32
            if (types >> t) & 1:
                rows.append(((hi << 64) | lo, us))
    lines = open(xcout).read().split("SOLUTION")[1].splitlines()[1:]
    sel = [int(l.split("row")[1]) for l in lines if l.strip().startswith("type")]
    comps = []
    for mask, us in (rows[r] for r in sel):
        cnt = {}
        for i in range(P):                   # inst.cpp encodes usage most significant first from (12, -), so the
            for d in (0, 1):                 # least significant digit is (point 0, D = +1)
                cnt[(i, d)] = us % 3
                us //= 3
        flags = [(i, d) for (i, d), c in cnt.items() for _ in range(c)]
        prs = [PAIRS[k] for k in range(78) if (mask >> k) & 1]
        comps.append((sorted({x for pr in prs for x in pr}), flags, prs))
    return comps


def hyper(u, v, c):
    return np.array([[int((u[0] * x[0] + u[1] * x[1] + v[0] * y[0] + v[1] * y[1] - c) % 3 == 0) for y in PTS9]
                     for x in PTS9], np.int64)


def t_dir(wd):
    return np.array([[int(any(((k * wd[0]) % 3, (k * wd[1]) % 3) == ((x[0] - y[0]) % 3, (x[1] - y[1]) % 3)
                              for k in range(3))) for y in PTS9] for x in PTS9], np.int64)


def build(comps):
    flag_class = {}
    for i in range(P):
        fl = [(ci, d) for ci, (pts, flags, prs) in enumerate(comps) for (pi, d) in flags if pi == i]
        plus, minus = [ci for ci, d in fl if d == 0], [ci for ci, d in fl if d == 1]
        assert len(plus) == 2 and len(minus) == 2
        flag_class.update({(plus[0], i): "A1", (plus[1], i): "A2", (minus[0], i): "B1", (minus[1], i): "B2"})
    rook = t_dir((1, 1)) + t_dir((1, 2)) - 2 * I9
    n = 9 * P
    C = np.zeros((n + 1, n + 1), np.int64)
    for i in range(P):
        C[9 * i:9 * i + 9, 9 * i:9 * i + 9] = J9 - I9 - 2 * rook
    for ci, (pts, flags, prs) in enumerate(comps):
        assert len(pts) == 4                 # K_4 components only
        Dv = {p: (1 if d == 0 else -1) for p, d in flags}
        E = list(itertools.combinations(pts, 2))
        tree = [(pts[0], q) for q in pts[1:]]
        non = [e for e in E if e not in tree]
        found = None
        for ph in itertools.product(range(3), repeat=len(non)):
            c = {e: 0 for e in tree}
            c.update(dict(zip(non, ph)))
            Phi = {}
            for e in E:
                v = chi(e[0] - e[1]) * W3 ** c[e]
                Phi[e], Phi[e[::-1]] = v, np.conj(v)
            W = np.array([[Dv[a] if a == b else 2 * Phi[(a, b)] for b in pts] for a in pts])
            if np.abs(W @ W - 13 * np.eye(4)).max() < 1e-9:
                found = Phi
                break
        assert found, f"no phases for component {ci}"
        for i, j in E:
            s = chi(i - j)
            cpp = int(round(np.angle(found[(i, j)] / s) / (2 * np.pi / 3))) % 3
            b = CLASSREP[flag_class[(ci, j)]]
            blk = -s * (J9 - 2 * hyper(CLASSREP[flag_class[(ci, i)]], ((-b[0]) % 3, (-b[1]) % 3), (-cpp) % 3))
            C[9 * i:9 * i + 9, 9 * j:9 * j + 9] = blk
            C[9 * j:9 * j + 9, 9 * i:9 * i + 9] = blk.T
    C[n, :n] = C[:n, n] = 1
    return C


def main():
    C = build(read_rows("insts.bin", sys.argv[1]))
    m = C.shape[0]
    assert np.array_equal(C, C.T) and not np.diag(C).any() and np.array_equal(C @ C, (m - 1) * np.eye(m, dtype=np.int64))
    perm = [m - 1] + list(range(m - 1))
    C = C[np.ix_(perm, perm)]
    s = np.concatenate(([1], np.where(C[0, 1:] < 0, -1, 1)))
    C = C * s[:, None] * s[None, :]
    text = "".join("".join({0: "0", 1: "+", -1: "-"}[int(x)] for x in r) + "\n" for r in C)
    name = f"conference118_{hashlib.sha256(text.encode()).hexdigest()[:16]}.txt"
    open(name, "w").write(text)
    print(f"C = C^T, zero diagonal, C^2 = {m - 1} I: checked; wrote {name}")


if __name__ == "__main__":
    main()

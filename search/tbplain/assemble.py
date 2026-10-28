#!/usr/bin/env python3
"""Core -> quaternary Hadamard matrix: exact check of the core, border (border.py), K, exact K K* = n I, certificate.

Takes four +-1 sequences A, B, C, D on G (a hit of n358/produce.py, a solution of n498/verify_sol.py, or an existing
certificate) and the structure ("358": G = Z2^4 x Z11, phi = (e2, e3), Q = Z2^2, T = the three nontrivial characters; "498":
G = Z2 x Z3 x Z41, phi = the Z3 coordinate, Q = Z3, T = all characters), and writes complex-hadamard/C<n>.json in the
format of the "szollosi-tbordered" certificates (verify/verify_complex.py): the sequences, the border (Xb, Rh, Ch) and the SHA-256
of the canonical text of K (n lines of n digits, digit e for the entry i^e).  Steps, all exact integer arithmetic:
(1) A and B symmetric and sum_k PAF_k(g) = -Gamma(phi(g)) over G by the group convolution (Gamma = 4 sum of the characters
of T, 4|G| at g = 0); (2) the compressions and the border; (3) K = [[Xb, R], [R', Kc]] with Kc = [[X, Y], [Y*, -X*]],
X = a X_A + conj(a) X_B, Y = a X_C + conj(a) X_D R, a = (1 + i)/2, R[a, (i, g)] = Rh[a, i q + phi(g)], R'[(i, g), b] =
Ch[i q + phi(g), b]; entries in {1, -1, i, -i}; K K* = n I (and K* K).  The statement is TBORDERED.md; verify_complex.py builds K
from the certificate independently.

    python3 assemble.py STRUCTURE CORE.json --out DIR [--index I] [--text]
      STRUCTURE  358 or 498 (102 and 18: the small analogues of the controls)
      CORE.json  {"seqs": [A, B, C, D]}, or a certificate ("sequences"), or a list of such records (--index picks one)
      --out DIR  writes DIR/complex-hadamard/C<n>.json (and C<n>.txt with --text); prints the SHA-256 of the matrix
"""
import argparse, hashlib, itertools, json, os, sys
import numpy as np
from border import Quot, border_search, finish

STRUCTURES = {
    "358": dict(mods=[2, 2, 2, 2, 11], quot=[2, 2], T=[[0, 1], [1, 0], [1, 1]], phi=lambda c: 2 * c[:, 2] + c[:, 3], aut=[]),
    "498": dict(mods=[2, 3, 41], quot=[3], T=[[0], [1], [2]], phi=lambda c: c[:, 1], aut=[[[1, 0, 0], [0, 1, 0], [0, 0, 10]]]),
    # the small analogues of the two structures (controls): N = 51 (order 102) and N = 9 (order 18)
    "102": dict(mods=[2, 2, 2, 2, 3], quot=[2, 2], T=[[0, 1], [1, 0], [1, 1]], phi=lambda c: 2 * c[:, 2] + c[:, 3], aut=[]),
    "18": dict(mods=[2, 3], quot=[3], T=[[0], [1], [2]], phi=lambda c: c[:, 1], aut=[]),
}


def group_developed(x, coords, mods):
    d = (coords[None, :, :] - coords[:, None, :]) % np.array(mods)
    idx = np.ravel_multi_index(tuple(np.moveaxis(d, 2, 0)), tuple(mods))
    return np.asarray(x, dtype=np.int64)[idx]          # X[g, h] = x(h - g)


def check_core(X, coords, mods, phi, Qt):
    """exact: A, B symmetric, sum_k PAF_k = target over G by the definition of the PAF"""
    m = len(coords)
    add = np.ravel_multi_index(tuple(((coords[:, None, :] + coords[None, :, :]) % np.array(mods)).transpose(2, 0, 1)), tuple(mods))
    neg = np.ravel_multi_index(tuple(((-coords) % np.array(mods)).T), tuple(mods))
    assert (X[0] == X[0][neg]).all() and (X[1] == X[1][neg]).all(), "A, B not symmetric"
    gamma = np.rint(4 * sum(Qt.chars[Qt.charkeys.index(tuple(t))] for t in Qt.T).real).astype(np.int64)      # Gamma(a)
    assert np.allclose(4 * sum(Qt.chars[Qt.charkeys.index(tuple(t))] for t in Qt.T).imag, 0)
    target = -gamma[phi]; target[0] = 4 * m
    tot = np.zeros(m, dtype=np.int64)
    for x in X:
        tot += np.array([int((x * x[add[g]]).sum()) for g in range(m)])
    assert (tot == target).all(), "sum of PAFs != target"


def assemble(structure, seq_strings, outdir=None, text=False):
    st = STRUCTURES[structure]; mods = st["mods"]
    m = int(np.prod(mods)); t = len(st["T"]); N = m + t
    coords = np.array(list(itertools.product(*[range(q) for q in mods])), dtype=np.int64)
    X = [np.array([1 if c == "+" else -1 for c in s], dtype=np.int64) for s in seq_strings]
    assert all(len(x) == m for x in X), "wrong length"
    Qt = Quot(st["quot"], st["T"]); q = Qt.q
    phi = np.asarray(st["phi"](coords), dtype=np.int64)
    check_core(X, coords, mods, phi, Qt)
    ys = [np.bincount(phi, weights=x, minlength=q).astype(np.int64) for x in X]
    r = border_search(Qt, ys)
    assert r is not None, "no border for this compression tuple"
    Rh, Ch = r; Xb = finish(Qt, ys, Rh, Ch)
    A, B, C, D = (group_developed(x, coords, mods) for x in X)
    neg = np.ravel_multi_index(tuple(((-coords) % np.array(mods)).T), tuple(mods))
    Rm = np.zeros((m, m), np.int64); Rm[np.arange(m), neg] = 1
    DR = D @ Rm
    Kcr = np.block([[(A + B) // 2, (C + DR) // 2], [(C.T + DR) // 2, -(A + B) // 2]])
    Kci = np.block([[(A - B) // 2, (C - DR) // 2], [-(C.T - DR) // 2, (A - B) // 2]])
    re, im = np.array([1, 0, -1, 0]), np.array([0, 1, 0, -1])
    cols = np.concatenate([phi, q + phi])
    Kr = np.block([[re[Xb], re[Rh][:, cols]], [re[Ch][cols, :], Kcr]])
    Ki = np.block([[im[Xb], im[Rh][:, cols]], [im[Ch][cols, :], Kci]])
    n = Kr.shape[0]; assert n == 2 * N
    assert np.all(np.abs(Kr) + np.abs(Ki) == 1), "entry not in {1, -1, i, -i}"
    I = np.eye(n, dtype=np.int64)
    assert np.array_equal(Kr @ Kr.T + Ki @ Ki.T, n * I) and not (Ki @ Kr.T - Kr @ Ki.T).any(), "K K* != n I"
    assert np.array_equal(Kr.T @ Kr + Ki.T @ Ki, n * I) and not (Kr.T @ Ki - Ki.T @ Kr).any(), "K* K != n I"
    e = np.where(Kr == 1, 0, np.where(Ki == 1, 1, np.where(Kr == -1, 2, 3)))
    txt = "".join("".join(map(str, row)) + "\n" for row in e)
    sha = hashlib.sha256(txt.encode()).hexdigest()
    cert = {"order": n, "construction": "szollosi-tbordered", "N": N, "group": mods, "quotient": st["quot"], "T": st["T"],
            "phi": [int(v) for v in phi], "automorphisms": st["aut"], "rows": [abs(int(x.sum())) for x in X],
            "sequences": ["".join("+" if v == 1 else "-" for v in x) for x in X],
            "border": {"Xb": Xb.tolist(), "Rh": Rh.tolist(), "Ch": Ch.tolist()}, "sha256": sha}
    if outdir:
        d = os.path.join(outdir, "complex-hadamard"); os.makedirs(d, exist_ok=True)
        json.dump(cert, open(os.path.join(d, "C%d.json" % n), "w"))
        if text: open(os.path.join(d, "C%d.txt" % n), "w").write(txt)
    return cert


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("structure", choices=sorted(STRUCTURES)); ap.add_argument("core"); ap.add_argument("--out", required=True)
    ap.add_argument("--index", type=int, default=0); ap.add_argument("--text", action="store_true")
    a = ap.parse_args()
    d = json.load(open(a.core))
    if isinstance(d, list): d = d[a.index]
    seqs = d.get("seqs") or d["sequences"]
    cert = assemble(a.structure, seqs, a.out, a.text)
    print("order %d: core exact, border found, K K* = K* K = %d I exact, entries in {1, -1, i, -i}" % (cert["order"], cert["order"]))
    print("sha256 of the canonical text of K: %s" % cert["sha256"])


if __name__ == "__main__":
    main()

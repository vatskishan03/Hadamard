#!/usr/bin/env python3
"""Lines (Walsh characters of E = Z2^4) of the H-invariant sequences of G = Z2^4 x Z11: label rows, line targets, line PAF.

Below the Z11 level the tower refines a sequence by its "lines".  The sixteen characters eta of E fall into six H-orbits
(one per row of the 6 x 6 matrix M[j][c] = sum over the E-class c of (-1)^(eta_j . e)); the first row is the weight row
w = (1, 3, 3, 3, 3, 3).  For a sequence x(c, z) the line j is the integer sequence a_j(z) = sum_c M[j][c] x(c, z) on Z11.
`line_matrix` returns the distinct rows of M (in the order used everywhere: rows[0] = weights, rows[1] = line 1, rows[2]
= line 2), `line_targets` the targets tau_j(s), s = 0..5, of the PAF of the lines (the Walsh transform over E of the target
of the sum of PAFs, restricted to the shifts s of Z11), and `paf_line` the PAF of integer sequences on Z11 for the shifts
s = 0..5.  Library module (used by engine.py and tower.c); running the file prints the rows and checks tau_0.

    python3 lines.py
"""
import itertools
import numpy as np


def line_matrix(S, RL):
    """distinct Walsh rows (one per character orbit): dict row -> list of eta; columns = E-classes of kind C (same labels as RL.cls['C'])"""
    cl = RL.cls["C"]; ncl = int(cl.max()) + 1
    ec = RL.ecoords
    rows = {}
    for eta in itertools.product(range(2), repeat=ec.shape[1]):
        chi = (-1) ** ((ec @ np.array(eta)) % 2)
        row = tuple(int(chi[cl == c].sum()) for c in range(ncl))
        rows.setdefault(row, []).append(eta)
    return rows


def line_targets(S, RL, etas):
    """tau[j][s] = sum over g = (e, s) of chi_eta_j(e) target(g), s = 0..5"""
    p = int(S.mods[-1]); ec = RL.ecoords
    G = S.coords
    e_of_g = RL.piE
    res = []
    for eta in etas:
        chi_e = (-1) ** ((ec @ np.array(eta)) % 2)
        t = []
        for s in range(p // 2 + 1):
            msk = G[:, -1] == s
            t.append(int((chi_e[e_of_g[msk]] * S.target[msk]).sum()))
        res.append(t)
    return np.array(res)                                     # (6, 6)


def paf_line(a):
    p = a.shape[1]
    return np.stack([(a * np.roll(a, -s, axis=1)).sum(1) for s in range(p // 2 + 1)], 1)


if __name__ == "__main__":
    from structure import Struct
    from rlevel import RLevel
    S = Struct("i"); RL = RLevel(S)
    rows = line_matrix(S, RL)
    print("distinct Walsh rows:", len(rows))
    for r, etas in rows.items(): print(r, "multiplicity", len(etas))
    tau = line_targets(S, RL, [v[0] for v in rows.values()])
    print("line targets tau_j(s):", tau.tolist())
    # tau_0 is the Z11 compression of the target; check against the direct sums over G
    tZ = [int(S.target[S.coords[:, -1] == s].sum()) for s in range(6)]
    assert tau[0].tolist() == tZ
    print("tau_0 == Z11-level target: OK")

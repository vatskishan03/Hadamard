"""The structure of the n = 955 core: row-sum patterns and exact compression tuples (the levels above the search of s1_*.py).

usage: python3 structure955.py
G = Z2^3 x Z29 (order 232), H = <(tau, u)> of order 7 (tau the Singer cycle of GF(8), u = 16 or 7 = 16^3, elements of order 7 of Z29*),
Q = Z2^3, T = the 7 nontrivial characters of Q (t = 7), N = 232 + 7 = 239,
core order 2N = 478 = (955 + 1)/2. Prints, for both u,
  1. the orbit counts and the target classes: sum_k PAF_k = 926 at 0, -28 on K0 minus {0} (K0 = 0 x Z29), +4 off K0;
  2. the row-sum patterns (a, b, c, d) with a^2 + b^2 + c^2 + d^2 = 954 = 4N - 2 allowed by the dynamic-programming counts:
     (7, 21, 8, 20) and (21, 21, 6, 6);
  3. for each pattern the exact join of the compression vectors over the quotient Q = Z2^3 (zq.feasible): one tuple each,
     (y(0), y(a)) of the four blocks. Pattern 2 = (21, 21, 6, 6): yA = yB = (0, 3), yC = yD = (-1, 1), the one used by the search
     (pattern 1 = (7, 21, 8, 20) has no tau-invariant border, see border8.py);
  4. the same over the quotient Z29: 3316 feasible tuples for pattern 2, the tuples of s1_tuples.py (CHECKED: same count).
"""
import time
from zq import *

tau = mat([[0, 0, 1], [1, 0, 1], [0, 1, 0]])
Tn = [k for k in itertools.product(range(2), repeat=3) if any(k)]
for uu in (16, 7):
    H = block_diag(tau, mat([[uu]]))
    st = ZS([2, 2, 2, 29], [H], phi=lambda c: c[:, :3] % 2, T=Tn, qmods=[2, 2, 2], name=f"Z2^3xZ29 u={uu}")
    print(st.info(), " t =", st.t, " target classes:", sorted(set(st.tvec.tolist())))
    pats = st.patterns(sumsq=954)
    print("  row-sum patterns (a, b, c, d) with sum of squares 954:", pats)
    L1 = np.zeros((3, 4), dtype=np.int64); L1[0, 0] = L1[1, 1] = L1[2, 2] = 1
    q1 = Quot(st, L1, [2, 2, 2])
    L2 = np.zeros((1, 4), dtype=np.int64); L2[0, 3] = 1
    q2 = Quot(st, L2, [29])
    for name, qt in (("Z2^3", q1), ("Z29", q2)):
        print("  quotient", name, ": fibre classes L", qt.info["L"][2], "A", qt.info["A"][2], " classes on Q", len(qt.qrep), " target on Q", qt.qTc.tolist())
        for p in pats:
            t0 = time.time()
            r = feasible(st, qt, p, verbose=False)
            if r is None:
                print("   ", p, "skipped (lists too large for the dense join)"); continue
            data, sols, wsum = r
            print("   ", p, "feasible compression tuples:", len(sols), " mass 10^%.2f" % (log10(wsum) if wsum > 0 else -1), " %.1fs" % (time.time() - t0))
            if name == "Z2^3":
                for s in sols[:3]:
                    print("        compressions (y(0), y(a) for a != 0) of A, B, C, D:", [data[k][0][s[k]].tolist() for k in range(4)])

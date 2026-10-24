"""Row-sum (trivial character) feasibility of an origin-defect shape, K = Z_m prime-like (non-origin orbits all of size h, p = (nu-1)/2 pairs).
 Necessary conditions (EXACT):  R^s_X = S_X + h * sum_o u_X(o)  ,  R^t_X = tau_X + h * sum_o v_X(o)   (sums of 2p units; for A: pair sums +-(1+i) resp. +-(1-i)):
 (R1) sum |R^s|^2 = q n/2 + (m-1) q z    (R2) sum |R^t|^2 = 4m - n/2 - (m-1) z    (R3) sum Im(R^s conj R^t) = 0    (R4) 2 sum Re(R^s conj eps) = sum|R^t|^2 - q n/2
 z in {0,1} (observed modes).  usage: python3 rowfeas.py m q h  -> for every origin cost 4-tuple of the identity, says feasible or not.
"""
import sys, itertools, math
from collections import defaultdict
U = [1, 1j, -1, -1j]
def block_options(kind, cost, q, A=False):
    """all (S, eps, tau, x0) for a block of type kind/cost at the origin; A: skew origin"""
    out = []
    for x0 in U:
        if not A and x0 != 1: continue
        if A and x0 not in (1, -1j): continue
        for x1 in U:
            for x2 in U:
                if A and abs(x2 - 1j * x1.conjugate()) > 1e-9: continue
                e = (x1 + x2) / 2; t = (x1 - x2) / 2
                c = (q - 1) * abs(e) ** 2 + 2 * (x0 * e.conjugate()).real
                if abs(c - cost) > 1e-9: continue
                ty = 'E' if abs(e) < 1e-9 else ('F' if abs(abs(e) - 1) < 1e-9 else 'M')
                out.append((x0 - e, e, t, ty))
    return out
def gaussian_vals(par, bound):
    """gaussian integers w=a+bi with a+b == par mod 2 and |w| <= bound"""
    r = []
    b = int(bound)
    for a in range(-b, b + 1):
        for c in range(-b, b + 1):
            if (a + c) % 2 == par % 2 and a * a + c * c <= bound * bound: r.append(complex(a, c))
    return r
def feasible(m, q, h, p, costs, z, verbose=False):
    """necessary conditions R1+R4 on s-row sums and R2 on t-row sums (cross condition R3 not used here)"""
    A = block_options('A', costs[0], q, A=True)
    opts = [A] + [block_options('F', c, q) for c in costs[1:]]
    if any(len(o) == 0 for o in opts): return (False, False)
    n2 = sum(abs(e) ** 2 for (S, e, t, ty) in [o[0] for o in opts])
    tgt1 = q * n2 + (m - 1) * q * z
    tgt2 = 4 * m - n2 - (m - 1) * z
    tgt4 = tgt2 - q * n2
    ws = gaussian_vals(0, 2 * p)
    Js = [J for J in range(-p, p + 1) if (J - p) % 2 == 0]
    # s-level: per block set of (4 ns, 4 jt)
    ssets = []; tsets = []
    for bi, o in enumerate(opts):
        s = set(); tt = set()
        for (S, e, t, ty) in o:
            Rs = [S + h * (1 + 1j) * J for J in Js] if bi == 0 else [S + h * w for w in ws]
            Rt = [t + h * (1 - 1j) * J for J in Js] if bi == 0 else [t + h * w for w in ws]
            for rs in Rs: s.add((int(round(4 * abs(rs) ** 2)), int(round(4 * 2 * (rs * e.conjugate()).real))))
            for rt in Rt: tt.add(int(round(4 * abs(rt) ** 2)))
        ssets.append(s); tsets.append(tt)
    # t-level
    reach = {0}
    for tt in tsets: reach = {a + b for a in reach for b in tt}
    t_ok = int(round(4 * tgt2)) in reach
    # s-level (norm, jtrace)
    reach = {(0, 0)}
    for ss in ssets: reach = {(a[0] + b[0], a[1] + b[1]) for a in reach for b in ss}
    s_ok = (int(round(4 * tgt1)), int(round(4 * tgt4))) in reach
    return (s_ok, t_ok)
if __name__ == '__main__':
    from oo_shapes import combos
    m, q, h = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    p = (m - 1) // h // 2
    for cb in combos(m, q):
        # cb sorted blocks B<=C<=D with A first ; try all orderings of the free three? the test is symmetric in B,C,D
        res = {z: feasible(m, q, h, p, cb, z) for z in (0, 1)}
        print(m, q, 'h', h, 'p', p, cb, res)

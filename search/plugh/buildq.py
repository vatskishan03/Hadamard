"""From a K-level solution (4 blocks (s,t,e) of a plain pk.Shape) build the G-level +-1 sequences A,B,C,D on G = K x F_q, check them directly
(entries +-1, A,B symmetric, sum PAF = 4|G| delta) and build Szollosi's quaternary matrix K = [[X, Y], [Y*, -X*]] of order 2|G| with
X = w X_A + w~ X_B, Y = w X_C + w~ X_D R,  w = (1+i)/2,  X_f[g,h] = f(h-g),  R[g,-g] = 1; check K K* = 2|G| I exactly (integers)."""
import numpy as np, sys, json
sys.path.insert(0, '.')
from idx import GF, Abel, factor_prime_power


def glevel_seq(sh, blk, F, G, q):
    """blk = (s,t,e,...) K-level arrays -> f on G (index k*q + a), a in F_q (GF index); classes by squares of F*"""
    s, t, e = [np.asarray(x, dtype=np.int64) for x in blk[:3]]
    x0 = s + e; x1 = e + t; x2 = e - t
    sq = {int(F.mul[x, x]) for x in range(1, q)}
    f = np.zeros(G.n, dtype=np.int64)
    for k in range(sh.N):
        for a in range(q):
            f[k * q + a] = x0[k] if a == 0 else (x1[k] if a in sq else x2[k])
    return f


def check_solution(sh, sol, verbose=True):
    q = sh.q
    p, ee = factor_prime_power(q)
    F = GF(p, ee)
    G = Abel(list(sh.Kmods) + [p] * ee)
    # note: sh.Kg element order is C order over Kmods, same as Abel's, so index k*q + a matches
    seqs = [glevel_seq(sh, b, F, G, q) for b in sol]
    n = G.n
    ok = all(set(np.unique(f)) <= {-1, 1} for f in seqs)
    symm = all((seqs[i][G.neg] == seqs[i]).all() for i in (0, 1))
    paf = sum(G.paf(f) for f in seqs)
    tgt = np.zeros(n, dtype=np.int64); tgt[0] = 4 * n
    pafok = bool((paf == tgt).all())
    # quaternary matrix
    idxm = np.arange(n)
    def Xf(f):
        # X_f[g,h] = f(h-g)
        hm = G.add[idxm[None, :], G.neg[idxm][:, None]]      # h - g
        return f[hm]
    XA, XB, XC, XD = [Xf(f) for f in seqs]
    R = np.zeros((n, n), dtype=np.int64); R[idxm, G.neg[idxm]] = 1
    # X = w XA + w~ XB : w = (1+i)/2  -> real = (A+B)/2, imag = (A-B)/2 ; w~ = (1-i)/2
    Xr = (XA + XB) // 2; Xi = (XA - XB) // 2
    DR = XD @ R
    Yr = (XC + DR) // 2; Yi = (XC - DR) // 2
    # K = [[X, Y],[Y*, -X*]]  ; Y* = conj transpose
    Kr = np.block([[Xr, Yr], [Yr.T, -Xr.T]])
    Ki = np.block([[Xi, Yi], [-Yi.T, Xi.T]])        # -X* imag part: -(conj X)^T imag = +Xi^T ; Y* imag = -Yi^T
    # K K* = (Kr + i Ki)(Kr^T - i Ki^T) = Kr Kr^T + Ki Ki^T + i (Ki Kr^T - Kr Ki^T)
    re = Kr @ Kr.T + Ki @ Ki.T; im = Ki @ Kr.T - Kr @ Ki.T
    entries_ok = bool(((np.abs(Kr) + np.abs(Ki)) == 1).all()) and bool((Kr * Ki == 0).all())
    had = bool((re == 2 * n * np.eye(2 * n, dtype=np.int64)).all() and (im == 0).all())
    if verbose:
        print(f'  G-level: entries +-1 {ok}, A/B symmetric {symm}, sum PAF = 4|G| delta {pafok}; quaternary K of order {2*n}: entries in {{+-1,+-i}} {entries_ok}, K K* = {2*n} I: {had}')
    return ok and symm and pafok and entries_ok and had, (seqs, (Kr, Ki))


if __name__ == '__main__':
    from pk import Shape
    from exact import exact_count
    Kmods = json.loads(sys.argv[1]); q = int(sys.argv[2]); hg = [(tuple(u) if isinstance(u, list) else u, ph) for u, ph in (json.loads(sys.argv[3]) if len(sys.argv) > 3 else [])]
    sh = Shape(Kmods, q, hg, mode='plain')
    tot, sols, info = exact_count(sh, want=3)
    print('K', Kmods, 'q', q, 'H', hg, 'K-level exact count', tot, 'sols', len(sols))
    for s in sols: print(check_solution(sh, s)[0])

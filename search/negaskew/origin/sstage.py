"""S-stage (decoupled mode z=0) existence/count for pattern (A skew M-type ; B,D F-type ; C E-type):
 Lock (per +-pair, 96 local states) + (S): sum_X P_{s_X}(k) = 0 for k != 0 (K = Z_m, H_K given).
 usage: python3 sstage.py m "[hK]" epsA xA epsB epsD   (eps given as complex strings e.g. '0.5+0.5j'), prints #Lock states, #(S)-solutions.
"""
import sys, os, itertools, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
from qshape import kshape
U = [1, 1j, -1, -1j]
def local_states(epsA, epsB, epsD):
    res = []
    for a in U:
        ua = (a, 1j * a.conjugate())
        for b1 in U:
            for b2 in U:
                for d1 in U:
                    for d2 in U:
                        tot = epsA * ua[0].conjugate() + epsA.conjugate() * ua[1] + epsB * b1.conjugate() + epsB.conjugate() * b2 + epsD * d1.conjugate() + epsD.conjugate() * d2
                        if abs(tot) < 1e-9: res.append((ua, (b1, b2), (d1, d2)))
    return res
def run(m, hK, epsA, xA, epsB, epsD, verbose=True, maxcomb=3e7):
    sh = kshape([m], hK); nu = sh['nu']; neg = sh['neg']; N = np.array(sh['N'])  # N[j][a][b]
    pairs = [o for o in range(1, nu) if o < neg[o]]
    p = len(pairs)
    LS = local_states(epsA, epsB, epsD)
    nL = len(LS)
    if verbose: print('nu', nu, 'pairs', p, 'local states', nL, 'combos', nL ** p)
    if nL ** p > maxcomb: return None
    # build u arrays for A,B,D over all combos
    idx = np.array(list(itertools.product(range(nL), repeat=p)))  # (ncomb, p)
    ncomb = len(idx)
    uA = np.zeros((ncomb, nu), complex); uB = np.zeros((ncomb, nu), complex); uD = np.zeros((ncomb, nu), complex)
    for t, o in enumerate(pairs):
        po = neg[o]
        a0 = np.array([LS[i][0][0] for i in range(nL)]); a1 = np.array([LS[i][0][1] for i in range(nL)])
        b0 = np.array([LS[i][1][0] for i in range(nL)]); b1 = np.array([LS[i][1][1] for i in range(nL)])
        d0 = np.array([LS[i][2][0] for i in range(nL)]); d1 = np.array([LS[i][2][1] for i in range(nL)])
        uA[:, o] = a0[idx[:, t]]; uA[:, po] = a1[idx[:, t]]
        uB[:, o] = b0[idx[:, t]]; uB[:, po] = b1[idx[:, t]]
        uD[:, o] = d0[idx[:, t]]; uD[:, po] = d1[idx[:, t]]
    # C: all unit assignments on orbits 1..nu-1
    allC = np.array(list(itertools.product(range(4), repeat=nu - 1)))
    Uarr = np.array(U)
    uC = np.zeros((len(allC), nu), complex); uC[:, 1:] = Uarr[allC]
    SA = xA - epsA; SB = 1 - epsB; SD = 1 - epsD; SC = 1.0
    # PAF at class lags j=1..nu-1: P_s(l_j) = sum_{a,b>=1} N[j][a][b] u_a conj(u_b) + S conj(u_{orb(l_j)}) + conj(S) u_{orb(-l_j)}
    def pafs(u, S):
        out = np.zeros((u.shape[0], nu - 1), complex)
        for j in range(1, nu):
            M = N[j][1:, 1:].astype(float)
            q = np.einsum('na,ab,nb->n', u[:, 1:], M, np.conj(u[:, 1:]))
            q = q + S * np.conj(u[:, j]) + np.conj(S) * u[:, neg[j]]
            out[:, j - 1] = q
        return out
    # careful: N[j][a][b] counts x in orbit a with x + l_j in orbit b : P(l)=sum_x u(x) conj u(x+l)
    PA = pafs(uA, SA); PB = pafs(uB, SB); PD = pafs(uD, SD); PC = pafs(uC, SC)
    tot = PA + PB + PD
    # exact integer keys: scale by 2
    def key(P):
        R = np.rint(2 * P.real).astype(np.int64); I = np.rint(2 * P.imag).astype(np.int64)
        assert np.max(np.abs(2 * P.real - R)) < 1e-6 and np.max(np.abs(2 * P.imag - I)) < 1e-6
        return R, I
    # need PC = -tot
    Rt, It = key(-tot); Rc, Ic = key(PC)
    from collections import defaultdict
    d = defaultdict(list)
    for i in range(len(allC)):
        d[(tuple(Rc[i]), tuple(Ic[i]))].append(i)
    nsol = 0; examples = []
    for i in range(ncomb):
        k = (tuple(Rt[i]), tuple(It[i]))
        if k in d:
            nsol += len(d[k])
            if len(examples) < 2: examples.append((i, d[k][0]))
    if verbose: print('(S)-solutions (A,B,D lock-consistent, C free):', nsol)
    return nL, p, nsol
if __name__ == '__main__':
    m = int(sys.argv[1]); hK = json.loads(sys.argv[2])
    c = lambda s: complex(s.replace(' ', ''))
    epsA = c(sys.argv[3]); xA = c(sys.argv[4]); epsB = c(sys.argv[5]); epsD = c(sys.argv[6])
    run(m, hK, epsA, xA, epsB, epsD)

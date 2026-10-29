"""Independent implementation (derived from the K-level algebra, not copied from rigid_s.py) of the rigid (A,M|E,E) K-level 'Lock + S' count:
 K = Z_m, H_K trivial.  Blocks X in {A,B,C,D}: s_X(k) = x0 (unit) off the origin for E-type, s_X(0) = x0 - eps_X;  A skew:  u_A(-k) = i conj u_A(k) (k != 0), x0_A in {1,-i}.
 A and B carry half defects at the origin (eps_X = (x1+x2)/2, cost (q-1)/2 + 2 Re(x0 conj eps)); C, D are E-type (eps = 0).
 Level-2 locking: Phi_A(k) + Phi_B(k) = 0 with Phi_X(k) = eps_X conj u_X(k) + conj eps_X u_X(-k)  (enumerated locally per +- pair, NOT assuming u_B = lambda u_A).
 Level-1 'S' equation:  sum_X P_{s_X}(k) = 0 for all k != 0,  P_s(k) = sum_j s(j) conj s(j+k).   C, D: all 4^(m-1) unit sequences with s(0) = 1 (global unit is a symmetry), no other reduction.
 usage: rigid_j.py m q [maxA] [variant]"""
import sys, itertools, time
import numpy as np
m = int(sys.argv[1]); q = int(sys.argv[2]); half = (m - 1) // 2
maxA = int(sys.argv[3]) if len(sys.argv) > 3 else 10**9
only = int(sys.argv[4]) if len(sys.argv) > 4 else None
U = [1, 1j, -1, -1j]
t0 = time.time()
rng = np.random.default_rng(12345)
R = rng.integers(1, 2**62, size=2 * half, dtype=np.int64) | 1
def cost(x0, x1, x2):
    e = (x1 + x2) / 2
    return (q - 1) * abs(e) ** 2 + 2 * (x0 * np.conj(e)).real
# origin variants
vA = []
for x0 in (1, -1j):
    for x1 in U:
        for x2 in U:
            if abs(x2 - 1j * np.conj(x1)) > 1e-9: continue      # skew at (0,f): u(0,-f) = i conj u(0,f)
            vA.append((x0, (x1 + x2) / 2, cost(x0, x1, x2)))
vB = []
for x1 in U:
    for x2 in U:
        # x0_B = 1 (global unit normalisation); defect types: all (x1,x2); require half defect or full defect, but exclude E (x2 = -x1: no defect)
        if abs(x2 + x1) < 1e-9: continue
        vB.append((1, (x1 + x2) / 2, cost(1, x1, x2)))
variants = []
seen = set()
for (x0A, eA, cA) in vA:
    for (x0B, eB, cB) in vB:
        if abs(cA + cB - 4 * m) < 1e-9:
            key = (round(x0A.real), round(x0A.imag), round(2 * eA.real), round(2 * eA.imag), round(2 * eB.real), round(2 * eB.imag))
            if key in seen: continue
            seen.add(key); variants.append((x0A, eA, x0B, eB))
print('m', m, 'q', q, 'variants (x0A, eA, eB):', [(v[0], v[1], v[3]) for v in variants], flush=True)
def paf_sc(S):       # S: (N, m) complex of Gaussian half-integers; return 4*P as integer array (N, half, 2)
    out = np.zeros((S.shape[0], half, 2), dtype=np.int64)
    for k in range(1, half + 1):
        P = (S * np.conj(np.roll(S, -k, axis=1))).sum(axis=1)
        out[:, k - 1, 0] = np.rint(4 * P.real); out[:, k - 1, 1] = np.rint(4 * P.imag)
    return out
def hkey(P):
    return (P.reshape(P.shape[0], -1) * R).sum(axis=1)
# C / D table: s(0)=1, others unit
idx = np.array(list(itertools.product(range(4), repeat=m - 1)), dtype=np.int8)
Cs = np.ones((idx.shape[0], m), complex); Cs[:, 1:] = np.array(U)[idx]
PC = paf_sc(Cs); kC = hkey(PC); order = np.argsort(kC); kCs = kC[order]
print('C table', Cs.shape[0], round(time.time() - t0, 1), flush=True)
total = 0
for vi, (x0A, eA, x0B, eB) in enumerate(variants):
    if only is not None and vi != only: continue
    # local states per +- pair
    local = []
    for a in U:
        a2 = 1j * np.conj(a)
        phiA = eA * np.conj(a) + np.conj(eA) * a2
        sol = [(b, b2) for b in U for b2 in U if abs(eB * np.conj(b) + np.conj(eB) * b2 + phiA) < 1e-9]
        local.append((a, a2, sol))
    nloc = sum(len(l[2]) for l in local)
    states = [(a, a2, b, b2) for (a, a2, sol) in local for (b, b2) in sol]
    print('variant', vi, 'x0A', x0A, 'eA', eA, 'eB', eB, 'local (a,b) states per pair:', nloc, flush=True)
    SA0 = x0A - eA; SB0 = x0B - eB
    cnt = 0; ncomb = len(states) ** half
    comb = itertools.product(range(len(states)), repeat=half)
    nA = 0
    for batch_start in range(0, min(ncomb, maxA), 4096):
        blk = list(itertools.islice(comb, 4096))
        if not blk: break
        N = len(blk)
        sA = np.zeros((N, m), complex); sB = np.zeros((N, m), complex)
        sA[:, 0] = SA0; sB[:, 0] = SB0
        for pos in range(half):
            for r, c in enumerate(blk):
                a, a2, b, b2 = states[c[pos]]
                k = pos + 1
                sA[r, k] = a; sA[r, m - k] = a2; sB[r, k] = b; sB[r, m - k] = b2
        T = -(paf_sc(sA) + paf_sc(sB))           # need PC + PD = T  (scaled by 4)
        kT = hkey(T)
        # D table = same as C table.  for each A-state: find d with kT - kD in kC
        for r in range(N):
            need = kT[r] - kC
            lo = np.searchsorted(kCs, need, 'left'); hi = np.searchsorted(kCs, need, 'right')
            nz = np.nonzero(hi - lo)[0]
            for d in nz:
                for ci in order[lo[d]:hi[d]]:
                    if np.array_equal(PC[ci] + PC[d], T[r]): cnt += 1
        nA += N
    print('  variant', vi, 'A-states', nA, 'of', ncomb, '(Lock+S) solutions:', cnt, ' t=', round(time.time() - t0, 1), flush=True)
    total += cnt
print('TOTAL', total)

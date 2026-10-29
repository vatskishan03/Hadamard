"""Numerical and exact checks of the statements proved in PROOFS.md (search/NOTES.md section 92): P1 to P7.  Every check tests a proved statement on real data or by exact integer computation; it does not replace the proof.
usage: python3 proofs_check.py [P1 P2 ...]        (default: all; needs ./tri2 built, see README.md; about 15 seconds)
  P1  parity of the PAF (Re odd, Im even) on random unit sequences over several odd moduli (also composite)
  P2  row sums |R|^2 = (4v-3, 1, 1, 1) on the exact solution sets of v = 29, 67 (tri2), on the recorded sample of v = 379 and the reduced list of v = 277; at v = 7 (k = 3) there are two patterns
  P3  the skew blocks: the bijection c -> (sigma, sigma'), the sizes of the lists by direct enumeration (v = 29 ... 631), the A of the solution sets lie in the list
  P4  conj U = (1+i) - iU at every t, |U|^2 = 2r^2 + 2r + 1, sum rho = +-(e - 2) or +-e, sum rho^2 = e(2v-7) or e(2v+1); the arithmetic of P4(c) for every k = 3 mod 4 up to 10^6
  P5  closure of the exact sets of v = 29, 67 under the group Z_e x tilde^3 x S_3 x negA, invariance of the pair sums; the 14 clusters of the 379 sample are the orbits of the group of order 32
  P6  the period form and the Gauss-sum form of the equations on a solution of v = 29, 67 and of v = 379
  P7  the lattice of the PAF vectors: even differences, rank e - 1, index 2^(e-1) in Z^(e-1) for v = 29, 67, 191, 277, 379, 631, 947, 1129 (exact integer elimination)"""
import sys, os, math, subprocess, itertools, collections, random
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from trilib import setup, expand, load_sols, I4
from tricore import Fam
import est_all as EA
from sympy import isprime

FAILS = []
def check(name, ok, info=''):
    print(('  ok   ' if ok else '  FAIL ') + name + (('   [' + info + ']') if info else ''), flush=True)
    if not ok: FAILS.append(name)

def tri2_sols(v, k, g, mode, M, tag):
    """all normalised solutions of (v, k) by tri2 (red=0), as in qskew's count: list of 4-tuples of blocks (u0 exponent, tuple of the e coset exponents)"""
    os.makedirs(os.path.join(HERE, 'runs'), exist_ok=True); fn = os.path.join(HERE, 'runs', 'pc_%s.out' % tag)
    with open(fn, 'w') as f:
        subprocess.run([os.path.join(HERE, 'tri2'), str(v), str(k), str(g), str(mode), '0', str(M), '0', '1', str(M), '0', '1000000', '1'], stdout=f, stderr=subprocess.DEVNULL, check=True)
    return load_sols(fn, v, k)[1]

def rowsum2(k, b):  return int(round(abs(I4[b[0]] + k * sum(I4[c] for c in b[1])) ** 2))
def paf_vec(F, b):
    U = np.zeros(F.e + 1, dtype=np.int64); U[0] = b[0]; U[1:] = b[1]
    return np.array(F.paf(I4[U][None, :])[0])
def shift(e, b, a): return (b[0], tuple(b[1][(j + a) % e] for j in range(e)))
def tilde(e, b):    return (b[0], tuple((-b[1][(j + e // 2) % e]) % 4 for j in range(e)))
def negA(b):        return ((3 if b[0] == 0 else 0), tuple((c + 2) % 4 for c in b[1]))
def image(e, s, a, na, tl, perm):
    A = shift(e, s[0], a); A = negA(A) if na else A
    fr = [shift(e, tilde(e, s[1 + i]) if tl[i] else s[1 + i], a) for i in range(3)]
    return (A,) + tuple(fr[p] for p in perm)

# ---------------------------------------------------------------------------------------------------------------- P1
def P1():
    print('P1  parity of the PAF')
    rng = np.random.default_rng(1); bad = 0; n = 0
    for v in (7, 15, 21, 29, 67, 101):
        for _ in range(300):
            u = I4[rng.integers(0, 4, size=v)]
            for s in range(1, v):
                p = np.sum(np.roll(u, -s) * np.conj(u)); n += 1
                if int(round(p.real)) % 2 != 1 or int(round(p.imag)) % 2 != 0: bad += 1
    check('P1: Re P(s) odd and Im P(s) even for every s != 0', bad == 0, '%d (sequence, shift) pairs over v = 7, 15, 21, 29, 67, 101' % n)

# ---------------------------------------------------------------------------------------------------------------- P2
def P2():
    print('P2  row sums')
    for (v, k, g, mode, M) in ((29, 7, 2, 7, 7), (67, 11, 2, 3, 31)):
        d = setup(v, k); sols = tri2_sols(v, k, g, mode, M, str(v)); pats = collections.Counter(tuple(rowsum2(k, b) for b in s) for s in sols)
        check('P2 at v = %d (k = %d, k mod 8 = %d): every solution has |R|^2 = (4v-3, 1, 1, 1) in some order' % (v, k, k % 8), len(sols) > 0 and all(sorted(p) == sorted((4 * v - 3, 1, 1, 1)) for p in pats), '%d solutions, patterns %s' % (len(sols), dict(pats)))
        if k % 8 == 3: check('   A is the big block (k = 3 mod 8)', all(rowsum2(k, s[0]) == 4 * v - 3 for s in sols))
        else:          check('   A is small and one free block is big (k = 7 mod 8)', all(rowsum2(k, s[0]) == 1 for s in sols))
    for (v, k, fn) in ((379, 27, 'data/s379.out'), (277, 23, 'data/c277_red.out')):
        d, sols = load_sols(os.path.join(HERE, fn), v, k); pats = collections.Counter(tuple(rowsum2(k, b) for b in s) for s in sols)
        check('P2 on the recorded %s' % fn, len(sols) > 0 and all(sorted(p) == sorted((4 * v - 3, 1, 1, 1)) for p in pats), '%d solutions, patterns %s' % (len(sols), dict(pats)))
    # v = 7 (k = 3): exhaustive over all 4-tuples
    v, k = 7, 3; F = Fam(v, k)
    allA = [(u0, (c0, (1 - c0) % 4)) for u0 in (0, 3) for c0 in range(4)]       # c_1 = i conj(c_0)
    free = [(0, (a, b)) for a in range(4) for b in range(4)]
    PA = [paf_vec(F, b) for b in allA]; PF = [paf_vec(F, b) for b in free]; pats = collections.Counter()
    for ia in range(len(allA)):
        for ib, ic, idd in itertools.product(range(len(free)), repeat=3):
            if not (PA[ia] + PF[ib] + PF[ic] + PF[idd]).any():
                pats[tuple(sorted((rowsum2(k, allA[ia]), rowsum2(k, free[ib]), rowsum2(k, free[ic]), rowsum2(k, free[idd])), reverse=True))] += 1
    check('P2 needs k >= 7: at v = 7 (k = 3) the patterns are (25,1,1,1) and (13,13,1,1)', set(pats) == {(25, 1, 1, 1), (13, 13, 1, 1)}, str(dict(pats)))

# ---------------------------------------------------------------------------------------------------------------- P3
def P3():
    print('P3  the skew blocks')
    sig = {}
    for x in range(4):
        c = I4[x]; s1 = (c + 1j * np.conj(c)) / (1 + 1j); s2 = (c - 1j * np.conj(c)) / (1 - 1j)
        assert abs(s1.imag) < 1e-12 and abs(s2.imag) < 1e-12
        sig[x] = (int(round(s1.real)), int(round(s2.real)))
    check("P3: c -> (sigma, sigma') is a bijection: 1 -> (+,+), i -> (+,-), -1 -> (-,-), -i -> (-,+)", sig == {0: (1, 1), 1: (1, -1), 2: (-1, -1), 3: (-1, 1)}, str(sig))
    check("P3: c = ((1+i) sigma + (1-i) sigma') / 2", all(abs(((1 + 1j) * p + (1 - 1j) * q) / 2 - I4[x]) < 1e-12 for x, (p, q) in sig.items()))
    for (v, k) in ((29, 7), (67, 11), (191, 19), (277, 23), (379, 27), (631, 35)):
        e = (v - 1) // k; h = e // 2; cnt = 0
        arr = np.array(list(itertools.product(range(4), repeat=h)), dtype=np.int64)        # c_0..c_{h-1}; c_{j+h} = i conj c_j has exponent 1 - x
        S = I4[arr].sum(1) + I4[(1 - arr) % 4].sum(1)
        for u0 in (0, 3):
            want = (1 + 1j) * (1 if u0 == 0 else -1) if k % 8 == 3 else 0
            cnt += int(np.sum(np.abs(S - want) < 1e-9))
        formula = 2 * math.comb(h, (h + 1) // 2) * 2 ** h if k % 8 == 3 else 2 * math.comb(h, h // 2) * 2 ** h
        check('P3 at v = %d (h = %d): |A| = %d by enumeration; |L| = C(e,e/2)^2 = %d, |Big| = C(e,e/2+1) C(e,e/2) = %d' % (v, h, cnt, math.comb(e, e // 2) ** 2, math.comb(e, e // 2 + 1) * math.comb(e, e // 2)), cnt == formula, 'formula %d' % formula)
    for (v, k, g, mode, M) in ((29, 7, 2, 7, 7), (67, 11, 2, 3, 31)):
        d = setup(v, k); e = d['e']; h = e // 2; sols = tri2_sols(v, k, g, mode, M, str(v))
        ok = all(b[0] in (0, 3) and all(b[1][(j + h) % e] == (1 - b[1][j]) % 4 for j in range(e)) for b in (s[0] for s in sols))
        check('P3: every A of the exact set of v = %d satisfies c_{j+h} = i conj(c_j), u0 in {1, -i}' % v, ok)

# ---------------------------------------------------------------------------------------------------------------- P4
def P4():
    print('P4  the spectrum of a skew block')
    for (v, k, g, mode, M) in ((29, 7, 2, 7, 7), (67, 11, 2, 3, 31)):
        d = setup(v, k); e = d['e']; sols = tri2_sols(v, k, g, mode, M, str(v)); As = sorted(set(s[0] for s in sols)); errs = 0.0; bad = 0; xs = np.arange(v)
        big = (k % 8 == 3); okS = True; okSS = True
        for b in As:
            u = expand(d, b); Uall = np.array([np.sum(u * np.exp(2j * np.pi * t * xs / v)) for t in range(v)])
            errs = max(errs, float(np.abs(np.conj(Uall) - ((1 + 1j) - 1j * Uall)).max()))
            r = (Uall - 1) / (1 + 1j); bad += int((np.abs(r.imag) > 1e-9).sum()); bad += int((np.abs(np.abs(Uall) ** 2 - (2 * r.real ** 2 + 2 * r.real + 1)) > 1e-8).sum())
            rho = np.array([(1 - 1j) * Uall[pow(g, m, v)] + 1j for m in range(e)])
            eps = 1 if b[0] == 0 else -1
            okS &= abs(rho.sum() - eps * (e - 2 if big else e)) < 1e-8
            okSS &= abs((rho ** 2).sum() - (e * (2 * v - 7) if big else e * (2 * v + 1))) < 1e-6 and np.abs(rho.imag).max() < 1e-9
        check('P4(a) at v = %d: conj U(t) = (1+i) - iU(t) at every t (max error %.1e); U = 1 + (1+i) r with r real; |U|^2 = 2r^2 + 2r + 1' % (v, errs), errs < 1e-9 and bad == 0, '%d distinct A' % len(As))
        check('P4(b) at v = %d: rho_m = (1-i) F(m) + i real, sum rho_m = %s e%s, sum rho_m^2 = %s' % (v, '+-(' if big else '+-', ' - 2)' if big else '', 'e(2v-7)' if big else 'e(2v+1)'), okS and okSS, 'e = %d' % e)
    bad = []
    for kk in range(7, 10 ** 6 + 1, 4):
        v = 1 + kk * (kk + 1) // 2; N = 2 * v - 7; assert N == kk * kk + kk - 5
        if math.isqrt(N) ** 2 == N: bad.append((kk, '2v-7 square'))
        if isprime(v) and N % v == 0 and math.isqrt(N // v) ** 2 == N // v: bad.append((kk, '2v-7 = v * square'))
        N1 = 2 * v + 1
        if math.isqrt(N1) ** 2 == N1: bad.append((kk, '2v+1 square'))
        if isprime(v) and N1 % v == 0 and math.isqrt(N1 // v) ** 2 == N1 // v: bad.append((kk, '2v+1 = v * square'))
    check('P4(c): for k = 3 mod 4, 7 <= k <= 10^6: 2v-7 = k^2+k-5 and 2v+1 are neither a square nor v times a square', not bad, str(bad[:3]))
    sq = [kk for kk in range(2, 5000) if math.isqrt(kk * kk + kk - 5) ** 2 == kk * kk + kk - 5]
    check('P4(c): k^2 + k - 5 is a square only for k = 2, 5 (identity (2k+1)^2 - 4n^2 = 21); at v = 7 (k = 3) 2v - 7 = 7 = v * 1^2 is the exception', sq == [2, 5] and 2 * 7 - 7 == 7, 'squares at k = %s' % sq)

# ---------------------------------------------------------------------------------------------------------------- P5
def P5():
    print('P5  symmetries and clusters')
    for (v, k, g, mode, M) in ((29, 7, 2, 7, 7), (67, 11, 2, 3, 31)):
        d = setup(v, k); e = d['e']; F = Fam(v, k); sols = tri2_sols(v, k, g, mode, M, str(v)); Sset = set(sols)
        def is_sol(s):
            return (not sum(paf_vec(F, b) for b in s).any()) and s[0][0] in (0, 3) and all(s[0][1][(j + e // 2) % e] == (1 - s[0][1][j]) % 4 for j in range(e)) and all(b[0] == 0 for b in s[1:])
        random.seed(5); smp = random.sample(sols, 40); bad = 0; nin = 0; nimg = 0
        for s in smp:
            for a in range(e):
                for na in (0, 1):
                    for tl in itertools.product((0, 1), repeat=3):
                        for perm in itertools.permutations(range(3)):
                            img = image(e, s, a, na, tl, perm); nimg += 1
                            if not is_sol(img): bad += 1
                            if img in Sset: nin += 1
        check('P5 at v = %d: the %d images of 40 solutions under the %d elements of Z_e x tilde^3 x S_3 x negA are solutions' % (v, nimg, 96 * e), bad == 0)
        if mode == 3: check('   (mode 3: the exact list is the whole normalised set, so every image lies in it)', nin == nimg, '%d of %d' % (nin, nimg))
        else: check('   (mode 7: the images with the big block at B and S = 1+i, %d of them, lie in the list)' % nin, nin > 0)
        okp = True
        for s in smp[:10]:
            L0 = paf_vec(F, s[0]) + paf_vec(F, s[1]); R0 = paf_vec(F, s[2]) + paf_vec(F, s[3])
            for na, tB, tC, tD, sw in itertools.product((0, 1), repeat=5):
                A = negA(s[0]) if na else s[0]; Bq = tilde(e, s[1]) if tB else s[1]; Cq = tilde(e, s[2]) if tC else s[2]; Dq = tilde(e, s[3]) if tD else s[3]
                if sw: Cq, Dq = Dq, Cq
                okp &= (paf_vec(F, A) + paf_vec(F, Bq) == L0).all() and (paf_vec(F, Cq) + paf_vec(F, Dq) == R0).all()
        check('   the left sum PAF(A)+PAF(B) and the right sum PAF(C)+PAF(D) are invariant under tilde_B, tilde_C, tilde_D, swap(C,D), negA (exact integer vectors)', bool(okp))
    d, sols = load_sols(os.path.join(HERE, 'data/s379.out'), 379, 27); e = d['e']
    rs = [l.split()[2] for l in open(os.path.join(HERE, 'data/s379.out')) if l.startswith('SOL')]
    byr = collections.defaultdict(list)
    for r, s in zip(rs, sols): byr[r].append(s)
    ok = True; sizes = set()
    for r, lst in byr.items():
        s = lst[0]; orb = set()
        for na, tB, tC, tD, sw in itertools.product((0, 1), repeat=5):
            A = negA(s[0]) if na else s[0]; fr = [tilde(e, s[1 + i]) if t else s[1 + i] for i, t in enumerate((tB, tC, tD))]
            if sw: fr[1], fr[2] = fr[2], fr[1]
            orb.add((A,) + tuple(fr))
        sizes.add(len(orb)); ok &= (orb == set(lst))
    check('P5: each of the %d hit partitions of the 379 sample is exactly one orbit of the group of order 32 (tilde B, C, D; swap C, D; A -> -A)' % len(byr), ok and sizes == {32} and len(sols) == 32 * len(byr), '%d solutions, orbit sizes %s' % (len(sols), sorted(sizes)))

# ---------------------------------------------------------------------------------------------------------------- P6
def P6():
    print('P6  period form and Gauss-sum form')
    cases = []
    for (v, k, g, mode, M) in ((29, 7, 2, 7, 7), (67, 11, 2, 3, 31)): cases.append((v, k, (lambda L: L[len(L) // 2])(tri2_sols(v, k, g, mode, M, str(v)))))
    cases.append((379, 27, load_sols(os.path.join(HERE, 'data/s379.out'), 379, 27)[1][0]))
    for (v, k, sol) in cases:
        d = setup(v, k); e = d['e']; h = e // 2; g = d['g']; F = Fam(v, k); eta = F.eta; zeta = np.exp(2j * np.pi / v); xs = np.arange(v)
        ok1 = np.allclose(eta[(np.arange(e) + h) % e], np.conj(eta)) and abs(eta.sum() + 1) < 1e-9
        Fd = []; Fp = []
        for b in sol:
            u = expand(d, b); Fd.append(np.array([np.sum(u * zeta ** ((pow(g, m, v) * xs) % v)) for m in range(e)]))
            Fp.append(np.array([I4[b[0]] + sum(I4[b[1][j]] * eta[(j + m) % e] for j in range(e)) for m in range(e)]))
        check('P6 at v = %d: eta_{j+h} = conj eta_j, sum_j eta_j = -1, and F(m) = u0 + sum_j c_j eta_{j+m} equals the direct spectrum' % v, ok1 and all(np.allclose(a, b) for a, b in zip(Fd, Fp)))
        tot = sum(np.abs(f) ** 2 for f in Fd); check('   sum_X |F_X(m)|^2 = 4v for every m', np.allclose(tot, 4 * v), 'max error %.1e' % np.abs(tot - 4 * v).max())
        X, Y = eta.real, eta.imag; sU = np.zeros(h); sIm = np.zeros(h); okp = True
        for b, f in zip(sol, Fd):
            c = I4[np.array(b[1])]; u0 = I4[b[0]]; dj = c[:h] + c[h:]; dpj = c[:h] - c[h:]
            Um = np.array([u0 + sum(dj[j] * X[(j + m) % e] for j in range(h)) for m in range(h)]); Vm = np.array([sum(dpj[j] * Y[(j + m) % e] for j in range(h)) for m in range(h)])
            okp &= np.allclose(f[:h], Um + 1j * Vm) and np.allclose(f[h:], Um - 1j * Vm)
            sU += np.abs(Um) ** 2 + np.abs(Vm) ** 2; sIm += (Um * np.conj(Vm)).imag
        check('   F(m) = U_m + iV_m, F(m+h) = U_m - iV_m; sum_X (|U|^2+|V|^2) = 4v and sum_X Im(U conj V) = 0 for m in Z_h', okp and np.allclose(sU, 4 * v) and np.allclose(sIm, 0), 'max errors %.1e, %.1e' % (np.abs(sU - 4 * v).max(), np.abs(sIm).max()))
        b = sol[0]; c = I4[np.array(b[1])]; u0 = I4[b[0]]
        sg = np.array([int(round(((c[j] + 1j * np.conj(c[j])) / (1 + 1j)).real)) for j in range(h)]); sp = np.array([int(round(((c[j] - 1j * np.conj(c[j])) / (1 - 1j)).real)) for j in range(h)])
        x = np.array([sum(sg[j] * X[(j + m) % e] for j in range(h)) for m in range(h)]); y = np.array([sum(sp[j] * Y[(j + m) % e] for j in range(h)) for m in range(h)])
        Um = np.array([u0 + sum((c[j] + c[j + h]) * X[(j + m) % e] for j in range(h)) for m in range(h)]); Vm = np.array([sum((c[j] - c[j + h]) * Y[(j + m) % e] for j in range(h)) for m in range(h)])
        okA = np.allclose(Um, u0 + (1 + 1j) * x) and np.allclose(Vm, (1 - 1j) * y) and (u0 != 1 or np.allclose((Um * np.conj(Vm)).imag, y * (2 * x + 1)))
        check('   skew block A: U_m = u0 + (1+i) x_m, V_m = (1-i) y_m, Im(U conj V) = y(2x+1) for u0 = 1', bool(okA))
        om = np.exp(2j * np.pi / e); cos = d['cos']
        G = np.array([sum(om ** ((l * cos[x_]) % e) * zeta ** x_ for x_ in range(1, v)) for l in range(e)])
        okG = abs(G[0] + 1) < 1e-9 and np.allclose(np.abs(G[1:]) ** 2, v)
        Bs = []
        for b in sol:
            c = I4[np.array(b[1])]; ch = np.array([c @ om ** (-l * np.arange(e)) for l in range(e)]) / e; Bl = ch * G; Bl[0] = I4[b[0]] - c.sum() / e; Bs.append(Bl)
        okF = all(np.allclose(f, np.array([sum(Bl[l] * om ** (-l * m) for l in range(e)) for m in range(e)])) for f, Bl in zip(Fd, Bs))
        paf_sum = sum(np.array([sum(Bl[(l + s) % e] * np.conj(Bl[l]) for l in range(e)) for s in range(e)]) for Bl in Bs)
        err = np.abs(paf_sum - 4 * v * (np.arange(e) == 0)).max()
        check('   Gauss-sum form: G_0 = -1, |G_l|^2 = v, F(m) = sum_l B(l) w^(-lm), sum_X PAF^{Z_e}_{B_X}(s) = 4v delta_s', okG and okF and err < 1e-6, 'max error %.1e' % err)

# ---------------------------------------------------------------------------------------------------------------- P7
def P7():
    print('P7  the lattice of the PAF vectors')
    for (v, k) in ((29, 7), (67, 11), (191, 19), (277, 23), (379, 27), (631, 35), (947, 43), (1129, 47)):
        F = Fam(v, k, seed=11); e, h = F.e, F.h; rng = F.rng; order = F.order; B = 800
        def s_askew(big):
            u0 = rng.integers(0, 2, size=B); npl = ((h + np.where(u0 == 0, 1, -1)) // 2) if big else np.full(B, h // 2)
            plus = np.take_along_axis(np.arange(h)[None, :] < npl[:, None], np.argsort(rng.random((B, h)), axis=1), axis=1)
            c = np.where(plus, rng.integers(0, 2, size=(B, h)), rng.integers(2, 4, size=(B, h)))
            ex = np.zeros((B, e + 1), dtype=np.int64); ex[:, 0] = np.where(u0 == 0, 0, 3); ex[:, 1:h + 1] = c; ex[:, h + 1:] = (1 - c) % 4; return I4[ex]
        def s_big1():
            a = rng.integers(0, h, size=B); b = h - 1 - a; idx = np.arange(e)[None, :]
            base = np.where(idx < (a + 1)[:, None], 0, np.where(idx < (2 * a + 1)[:, None], 2, np.where(idx < (2 * a + 1 + b + 1)[:, None], 1, 3)))
            ex = np.zeros((B, e + 1), dtype=np.int64); ex[:, 1:] = np.take_along_axis(base, np.argsort(rng.random(base.shape), axis=1), axis=1); return I4[ex]
        types = [s_askew(k % 8 == 3), F.s_small(B)] + ([] if k % 8 == 3 else [s_big1()])
        diffs = []
        for U in types:
            P = F.paf(U)[:, order]
            diffs += [(P[i] - P[0]) for i in range(1, B)]
        par = all(int(x) % 2 == 0 for dd in diffs for x in dd)
        basis = EA.int_basis(diffs); rk = len(basis); idx = 1
        for b in basis: idx *= abs(next(t for t in b if t))
        check('P7 at v = %d (e = %d): differences even; rank %d, index %d = 2^%d in Z^%d' % (v, e, rk, idx, e - 1, e - 1), par and rk == e - 1 and idx == 2 ** (e - 1))

if __name__ == '__main__':
    todo = sys.argv[1:] or ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7']
    for t in todo: globals()[t]()
    print('ALL CHECKS PASSED' if not FAILS else 'FAILED: ' + ', '.join(FAILS))
    sys.exit(1 if FAILS else 0)

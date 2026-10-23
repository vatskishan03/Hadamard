"""verify_all.py: all numerical checks of the octic-ring derivations (NOTES.md section 89), one PASS line each, ALL PASS at the end.  CPU: single thread,
about 2 minutes.  Ring isomorphism Z[x]/(x^{4v}+1) = Z[w][Z_v] (v = 1 .. 15: product, adjoint, row sum, det 1, explicit inverse), the 16-letter alphabet
(norm 8, moduli 4 +- 2 sqrt2), sum U U* = image of the negaperiodic autocorrelations with the peak 4n + sqrt2 NPAF(v), skew <=> palindromic, the row sums of the
exhaustive v = 3 solutions, the multiplier classification, and the theta-paired solutions.  Needs ./exsmall (cc -O2 -o exsmall exsmall.c).
usage: python3 verify_all.py"""
import numpy as np, random, itertools, subprocess, sys, json, math, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from octic import *
from dbl import mono_images, theta, adj
from mkcert import bits_to_seq, mkcert
random.seed(7); np.random.seed(7)
s2 = 2 ** .5
ok_all = True
def check(name, cond):
    global ok_all
    ok_all &= bool(cond); print(('PASS ' if cond else 'FAIL ') + name, flush=True)
# 1 ring isomorphism
for v in (1, 3, 5, 7, 9, 11, 13, 15):
    n = 4 * v; good = True
    for t in range(4):
        a = np.random.randint(-3, 4, n); b = np.random.randint(-3, 4, n)
        good &= np.array_equal(phi(nega_mul(a, b), v), gmul(phi(a, v), phi(b, v)))
        good &= np.array_equal(phi(nega_adj(a), v), gadj(phi(a, v)))
        good &= np.array_equal(zeta_eval(a, v), phi(a, v).sum(axis=0))
        good &= np.array_equal(psi(phi(a, v), v), a)
    M = np.array([phi(np.eye(n, dtype=np.int64)[k], v).reshape(-1) for k in range(n)])
    good &= round(np.linalg.det(M.astype(float))) == 1
    check(f'v={v}: phi ring hom (product, adjoint, DC), det=1, psi inverse', good)
# 2 alphabet
A = alphabet(); S = {tuple(x) for x in A}
vals = set()
for v in (3, 5, 7, 9, 11):
    for t in range(500):
        a = np.random.choice([-1, 1], 4 * v); U = phi(a, v)
        for g in range(v): vals.add(tuple(U[g]))
check('alphabet: all values of phi(+-1 sequence) lie in A, all 16 occur', vals == S)
check('alphabet closed under w, conj, -conj', all(tuple(zmul(c, zpow_w(1))) in S and tuple(zconj(c)) in S and tuple(-zconj(c)) in S for c in A))
mods = [abs(zcomplex(c)) ** 2 for c in A]
check('|c|^2 in {4+2sqrt2 (8 elements), 4-2sqrt2 (8 elements)}', sum(abs(m - (4 + 2 * s2)) < 1e-9 for m in mods) == 8 and sum(abs(m - (4 - 2 * s2)) < 1e-9 for m in mods) == 8)
def fnorm(c): return np.prod([sum(c[j] * np.exp(1j * np.pi * k * j / 4) for j in range(4)) for k in (1, 3, 5, 7)])
check('field norm of every alphabet element = 8 = N((1-w)^3)', all(abs(fnorm(c) - 8) < 1e-9 for c in A))
# 3 equation <-> NPAF, peak formula
for v in (3, 5, 7):
    n = 4 * v; good = True
    for t in range(5):
        seqs = [np.random.choice([-1, 1], n) for _ in range(4)]
        tot = sum(gmul(phi(s, v), gadj(phi(s, v))) for s in seqs)
        npaf = [sum(nega_paf(s, k) for s in seqs) for k in range(n)]
        # image of sum_s NPAF(s) x^s under phi
        img = phi(np.array(npaf, dtype=np.int64), v)
        good &= np.array_equal(tot, img)
        good &= np.array_equal(tot[0], [4 * n, 0, 0, 0] + np.array([0, 1, 0, -1]) * npaf[v]) if False else np.array_equal(tot[0], np.array([4 * n, npaf[v], 0, -npaf[v]]))
    check(f'v={v}: sum U U* = phi(sum NPAF x^s); peak = 4n + sqrt2 NPAF(v) (coordinates (4n, NPAF(v), 0, -NPAF(v)))', good)
# 4 skew <=> palindromic (exhaustive v=3, all 2^11 with a0 = 1 -> check both directions on all 2^12/2)
v = 3; n = 12; cnt = 0; agree = True
for bits in range(1 << (n - 1)):
    a = np.array([1] + [(-1) ** ((bits >> k) & 1) for k in range(n - 1)])
    pal = all(a[j] == a[n - j] for j in range(1, n))
    U = phi(a, v); T = U + gadj(U); sk = np.array_equal(T[0], [2, 0, 0, 0]) and not T[1:].any()
    agree &= (pal == sk); cnt += pal
check(f'v=3: palindromic <=> U + U* = 2 (all {1 << (n - 1)} sequences with a0=1), #palindromic = {cnt} = 2^{2 * v}', agree and cnt == 2 ** (2 * v))
# 5 DC on exhaustive v = 3 solutions (the SOL lines of ./exsmall 12)
os.makedirs('runs', exist_ok=True)
if not os.path.exists('runs/ex12.out'):
    with open('runs/ex12.out', 'w') as fo: subprocess.run(['./exsmall', '12', '1000000', '0'], stdout=fo, check=True)
sols = [list(map(int, l.split()[1:])) for l in open('runs/ex12.out') if l.startswith('SOL')]
good = True; pats = {}
for sol in sols[::50]:
    Rs = [zeta_eval(np.array(bits_to_seq(b, 12)), 3) for b in sol]
    P = sum(int((R ** 2).sum()) for R in Rs); Q = sum(int(R[0] * R[1] + R[1] * R[2] + R[2] * R[3] - R[0] * R[3]) for R in Rs)
    good &= (P == 48 and Q == 0 and Rs[0][0] == 1 and Rs[0][1] == Rs[0][3] and all(int(x) % 2 for R in Rs for x in R))
    key = tuple(sorted(tuple(int(x) for x in R) for R in Rs[1:]))
    pats[(tuple(int(x) for x in Rs[0]))] = pats.get(tuple(int(x) for x in Rs[0]), 0) + 1
check(f'v=3: on {len(sols[::50])} sampled exhaustive solutions sum P = 16v = 48, sum Q = 0, R_A = (1,r1,r2,r1), all coordinates odd; R_A patterns {sorted(pats)}', good)
# 6 compatible multipliers (re-run classification for a few v)
import mult_classify as t3
for v in (5, 7, 9, 15):
    N = 8 * v; seen = set(); has_tw = False
    for t in range(1, N, 2):
        if math.gcd(t, N) != 1: continue
        for s in range(N):
            for eps in (1, -1):
                if t == 1 and s == 0 and eps == 1: continue
                if t == 1 and eps == -1 and s == 4 * v: continue  # identity as -x^{4v} = 1
                c = t3.count_inv(v, [(eps, s, t)])
                if c >= 16:
                    seen.add((t % 8, t % v, s % v if False else 0))
    # every element admitting >= 16 invariant skew sequences must have t = 1 mod 8
    check(f'v={v}: elements (eps,s,t) with >=16 invariant skew sequences all have t = 1 (mod 8): {sorted(set(k[0] for k in seen))}', all(k[0] == 1 for k in seen))
# 7 Galois fixed points
def sigma(c, t):
    out = np.zeros(4, dtype=np.int64)
    for j in range(4): out += c[j] * zpow_w(j * t)
    return out
check('no alphabet element fixed by sigma_3, sigma_5, sigma_7 (nor w^k-twisted for sigma_3, sigma_5)', all(not np.array_equal(sigma(c, t), c) for c in A for t in (3, 5, 7)) and all(not np.array_equal(sigma(c, t), zmul(c, zpow_w(k))) for c in A for t in (3, 5) for k in range(8)))
# 8 doubling characterisation on v = 3
def even_odd(Y): return Y[0::2], Y[1::2]
def nega_paf_sum(seqs, n): return [sum(nega_paf(s, k) for s in seqs) for k in range(n)]
v = 3; n = 12; tested = 0; good = True; centre_pal = True
for sol in sols:
    S4 = [np.array(bits_to_seq(b, n)) for b in sol]
    A_, B_, C_, D_ = S4
    # find pairing
    found = None
    for Ap in (A_, adj(A_)):
        tA = theta(Ap)
        for X, (Y, Z) in ((B_, (C_, D_)), (C_, (B_, D_)), (D_, (B_, C_))):
            if not any(np.array_equal(im, X) for im in mono_images(tA)): continue
            for Yp, Zq in ((Y, Z), (Z, Y)):
                for Ypp in (Yp, adj(Yp)):
                    tY = theta(Ypp)
                    if any(np.array_equal(im, Zq) for im in mono_images(tY)): found = (Ap, Ypp); break
                if found: break
            if found: break
        if found: break
    if not found: continue
    Ap, Yp = found; a, b = even_odd(Ap); c, d = even_odd(Yp)
    tested += 1
    f1 = nega_paf_sum([a, b, c, d], n // 2)
    good &= (f1[0] == 4 * (n // 2) and all(x == 0 for x in f1[1:]))
    good &= (a[0] == 1 and all(a[j] == a[n // 2 - j] for j in range(1, n // 2)))
    centre_pal &= all(b[k] == b[n // 2 - 1 - k] for k in range(n // 2))
    if tested >= 400: break
check(f'v=3: {tested} theta-paired octic solutions give F1 quadruples (a skew, b centre-palindromic, sum NPAF_{{2v}} = 4*2v delta)', good and centre_pal and tested > 0)
print('ALL PASS' if ok_all else 'SOME FAIL')

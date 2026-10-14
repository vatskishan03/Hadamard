#!/usr/bin/env python3
"""Gaussian model for T3-twisted T-matrices of prime order t.

K = <H, mu> <= Z_t^*, [K:H] = 3, |H| = h.  Blocks (T1, T2, T3, T4) with T1 K-invariant (contains 0)
and T_{2+j}(x) = T2(mu^{-j} x).  Each K-orbit Q is either in T1 (sign s_Q) or of type T2@a
(a in {0,1,2}: T2 is supported on the H-coset g^{j+e*a}H of Q, sign s_Q).
Class values V(s) = sum_x [lab(x) = lab(x+s)] T(x) T(x+s) on <K,-1>-classes."""
import sys, math
import numpy as np
from sympy import primitive_root, binomial

def setup(t, h):
    g = primitive_root(t)
    d = 3 * h
    e = (t - 1) // d          # index of K = number of K-orbits; K = <g^e>, H = <g^(3e)>, mu = g^e
    q = (t - 1) // d
    # K-orbit index Q(x) and sub-index b(x) in {0,1,2} for x = g^k, k = j + e*b' ... use discrete log
    dlog = np.zeros(t, dtype=np.int64)
    x = 1
    for k in range(t - 1):
        dlog[x] = k; x = x * g % t
    Qof = np.full(t, -1); bof = np.full(t, -1)
    for x in range(1, t):
        k = dlog[x]
        Qof[x] = k % e                     # K-coset: k mod e
        bof[x] = (k // e) % 3              # which H-coset inside the K-coset
    # classes: <K,-1> cosets.  -1 = g^((t-1)/2).  <g^e, g^((t-1)/2)> = <g^gcd(e,(t-1)/2)>
    f = math.gcd(e, (t - 1) // 2)
    reps = [pow(g, c, t) for c in range(f)]
    return dict(t=t, h=h, g=g, e=e, Qof=Qof, bof=bof, reps=reps, C=f)

def class_values(S, typ, sgn):
    """typ[Q] in {-1 (T1), 0,1,2 (T2@a)}; sgn[Q] in +-1.  Batched: typ, sgn shape (B, e)."""
    t = S['t']; Qof, bof = S['Qof'], S['bof']
    B = typ.shape[0]
    lab = np.zeros((B, t), dtype=np.int64); val = np.zeros((B, t), dtype=np.int64)
    lab[:, 0] = 1; val[:, 0] = 1
    xs = np.arange(1, t)
    ty = typ[:, Qof[xs]]                    # (B, t-1)
    lab[:, 1:] = np.where(ty < 0, 1, 2 + (bof[xs][None, :] - ty) % 3)
    val[:, 1:] = sgn[:, Qof[xs]]
    out = np.zeros((B, len(S['reps'])), dtype=np.int64)
    x = np.arange(t)
    for c, s in enumerate(S['reps']):
        y = (x + s) % t
        same = lab[:, x] == lab[:, y]
        out[:, c] = (same * val[:, x] * val[:, y]).sum(1)
    return out

def rowsums(typ, sgn, h):
    r1 = 1 + 3 * h * np.where(typ < 0, sgn, 0).sum(1)
    r2 = h * np.where(typ >= 0, sgn, 0).sum(1)
    return r1, r2

def sample(S, k1, B, rng, r1t, r2t):
    """Random configurations with exactly k1 T1-orbits and the target row sums (rejection on signs)."""
    e = S['e']; h = S['h']
    typ = np.zeros((B, e), dtype=np.int64); sgn = np.zeros((B, e), dtype=np.int64)
    for b in range(B):
        perm = rng.permutation(e)
        T1 = perm[:k1]; T2 = perm[k1:]
        typ[b, T1] = -1; typ[b, T2] = rng.integers(0, 3, size=len(T2))
        # signs with prescribed sums: T1 signs sum (r1t-1)/(3h), T2 signs sum r2t/h
        s1 = (r1t - 1) // (3 * h); s2 = r2t // h
        p1 = (k1 + s1) // 2; p2 = (len(T2) + s2) // 2
        a = np.array([1] * p1 + [-1] * (k1 - p1)); rng.shuffle(a)
        c = np.array([1] * p2 + [-1] * (len(T2) - p2)); rng.shuffle(c)
        sgn[b, T1] = a; sgn[b, T2] = c
    return typ, sgn

if __name__ == '__main__':
    t, h, r1t, r2t = map(int, sys.argv[1:5])
    B = int(sys.argv[5]) if len(sys.argv) > 5 else 2000
    S = setup(t, h)
    e = S['e']; C = S['C']
    print(f"t={t} h={h} |K|={3*h} K-orbits={e} classes={C} g={S['g']}")
    rng = np.random.default_rng(1)
    s1 = (r1t - 1) // (3 * h); s2 = r2t // h
    totE = 0.0
    for k1 in range(0, e + 1):
        k2 = e - k1
        if (k1 + s1) % 2 or (k2 + s2) % 2 or abs(s1) > k1 or abs(s2) > k2: continue
        ncfg = float(binomial(e, k1) * binomial(k1, (k1 + s1) // 2) * binomial(k2, (k2 + s2) // 2) * 3 ** k2)
        typ, sgn = sample(S, k1, B, rng, r1t, r2t)
        V = class_values(S, typ, sgn)
        mu = V.mean(0); var = V.var(0)
        # lattice: gcd of differences per class
        lp = 0.0; terms = []
        for c in range(C):
            u = np.unique(V[:, c]); g = int(np.gcd.reduce(np.diff(u))) if len(u) > 1 else 0
            if g and (0 - u[0]) % g: lp = -math.inf; break
            terms.append(-0.5 * math.log10(2 * math.pi * var[c]) - mu[c] ** 2 / (2 * var[c]) / math.log(10) + math.log10(max(g, 1)))
        if lp == -math.inf:
            print(f"k1={k1}: lattice excludes 0"); continue
        terms.sort()
        lp = sum(terms[:-1])      # drop one (dependent) class
        lE = math.log10(ncfg) + lp
        totE += 10 ** lE
        print(f"k1={k1:2d} log10 cfg={math.log10(ncfg):6.2f} sd={np.sqrt(var).round(1).tolist()[:5]}.. lattice g={g} log10 P={lp:7.2f} log10 E={lE:6.2f}", flush=True)
    print(f"total log10 E = {math.log10(totE) if totE > 0 else -math.inf:.2f}")

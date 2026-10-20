#!/usr/bin/env python3
"""Negaperiodic Golay pairs of length v (NAF_p + NAF_q = 2v delta) invariant under a group H of multipliers
y -> y^k of Z[y]/(y^v + 1) (k odd, gcd(k, 2v) = 1): each sequence is constant up to sign on the signed orbits of H
(a structure is skipped when an orbit is forced to 0), the power-spectrum test |p(zeta)|^2 <= 2v at the odd 2v-th roots
of unity filters the list, and the two lists are joined exactly on the NAF (hadamard/MIYAMOTO_BLOCKS.md: an NG pair of
length (q - 1)/2 gives Miyamoto-block Hadamard matrices of order 4qt). Multiplier groups are tried by increasing number
of orbits, up to MAXO.

    python3 search/ng_mult.py v MAXO [MINO]      # e.g. 134 27 (QW(67), FRONTIER 1883); 14 8 finds a Paley pair
Hits go to $HADAMARD_WORK/ng_mult_hits.txt.
"""
import itertools, math, sys, time
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import workdir as W

def closure(gens, m):
    S = {1}; todo = [1]
    while todo:
        x = todo.pop()
        for g in gens:
            y = x * g % m
            if y not in S: S.add(y); todo.append(y)
    return sorted(S)

def signed_orbits(v, H):
    """orbits of H on positions 0..v-1 with signs: p[pos] = eps_orbit * sign. None if some orbit is forced to 0."""
    seen = {}
    orbits = []
    for j in range(v):
        if j in seen: continue
        orb = {j: 1}; todo = [j]; ok = True
        while todo:
            x = todo.pop()
            for k in H:
                e = (k * x) % (2 * v)
                y, s = (e, 1) if e < v else (e - v, -1)
                s *= orb[x]
                if y in orb:
                    if orb[y] != s: ok = False
                else:
                    orb[y] = s; todo.append(y)
        if not ok: return None
        for y in orb: seen[y] = len(orbits)
        orbits.append(orb)
    return orbits

def naf_fft(P, v):
    w = np.exp(1j * np.pi * np.arange(v) / v)
    F = np.fft.fft(P * w[None, :], axis=1)
    R = np.fft.ifft(np.abs(F) ** 2, axis=1)            # R[s] = sum_j conj(q_j) q_{j+s}... (convention checked below)
    return F, R

def search(v, H, log=print, maxlist=1 << 28):
    orbits = signed_orbits(v, H)
    if orbits is None: return None, 'forced zero'
    O = len(orbits)
    if 2 ** O > maxlist: return None, f'{O} orbits: too many'
    M = np.zeros((O, v))
    for i, orb in enumerate(orbits):
        for j, s in orb.items(): M[i, j] = s
    w = np.exp(1j * np.pi * np.arange(v) / v)
    Mw = M * w[None, :]
    Fm = np.fft.fft(Mw, axis=1)                          # F(eps) = eps @ Fm (linear)
    keepP, keepF = [], []
    nall = 2 ** (O - 1)
    CH = 1 << 15
    bits = np.arange(O)
    for c0 in range(0, nall, CH):
        idx = np.arange(c0, min(nall, c0 + CH), dtype=np.int64)
        E = 1.0 - 2.0 * ((idx[:, None] >> bits[None, :]) & 1)
        F = E @ Fm
        ok = np.all(np.abs(F) ** 2 <= 2 * v + 1e-6, axis=1)
        if ok.any():
            keepP.append(E[ok] @ M); keepF.append(F[ok])
    if not keepP: return None, f'{O} orbits, list 0 of {nall}'
    P = np.concatenate(keepP); F = np.concatenate(keepF)
    # NAF(s) = Re(conj(w^s) * sum_j q_j conj(q_{j+s})) via ifft of |F|^2
    R = np.fft.ifft(np.abs(F) ** 2, axis=1)            # R[s] = sum_j q_{j+s} conj(q_j)
    naf = np.rint((R * np.conj(w)[None, :]).real).astype(np.int64)
    # NAF values needed: s = 1..v-1 up to the H-action; use all s for safety but compress by unique columns
    cols = list(range(1, v))
    V = naf[:, cols]
    # exact join: V_p + V_q = 0
    rng = np.random.default_rng(1)
    wts = rng.integers(1, 2 ** 61, size=V.shape[1], dtype=np.int64)
    kp = (V * wts).sum(1); kq = (-V * wts).sum(1)
    order = np.argsort(kq); kqs = kq[order]
    lo = np.searchsorted(kqs, kp, 'left'); hi = np.searchsorted(kqs, kp, 'right')
    for i in np.flatnonzero(hi > lo):
        for j in order[lo[i]:hi[i]]:
            if np.array_equal(V[i] + V[j], np.zeros(V.shape[1], np.int64)):
                return (P[i].astype(np.int64), P[j].astype(np.int64)), f'{O} orbits, list {len(P)}'
    return None, f'{O} orbits, list {len(P)} of {nall}'

def is_ng(a, b):
    v = len(a)
    naf = lambda x: np.array([int(x[:v - s] @ x[s:]) - int(x[v - s:] @ x[:s]) for s in range(v)])
    t = naf(a) + naf(b)
    return t[0] == 2 * v and not t[1:].any()

def subgroups(v, maxorder=None):
    m = 2 * v
    units = [k for k in range(1, m) if math.gcd(k, m) == 1]
    out = {}
    for g in units:
        H = tuple(closure([g], m))
        out.setdefault(H, (g,))
    for (H1, g1), (H2, g2) in itertools.combinations(list(out.items()), 2):
        H = tuple(closure([g1[0], g2[0]], m))
        out.setdefault(H, (g1[0], g2[0]))
    return out

if __name__ == '__main__':
    v = int(sys.argv[1]); maxO = int(sys.argv[2]); minO = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    subs = subgroups(v)
    cand = []
    for H, g in subs.items():
        orb = signed_orbits(v, H)
        if orb is None: continue
        cand.append((len(orb), len(H), g, H))
    cand.sort()
    print(f"v = {v}: {len(cand)} admissible multiplier groups; orbit counts {sorted(set(c[0] for c in cand))[:30]}", flush=True)
    for O, h, g, H in cand:
        if O > maxO or O < minO: continue
        t = time.time()
        sol, info = search(v, H)
        print(f"  |H| = {h} gens {g}: {info}: {'FOUND' if sol else 'none'} [{time.time() - t:.1f}s]", flush=True)
        if sol:
            print('  NG check:', is_ng(*sol), flush=True)
            with open(W.hits('ng_mult_hits.txt'), 'a') as f:
                f.write(f"v={v} H={g}\n" + ''.join(''.join('+' if z > 0 else '-' for z in x) + '\n' for x in sol))
            break

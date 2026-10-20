"""Exhaustive Legendre-pair search over small multiplier structures.

A Legendre pair (A, B) of odd length l: +-1 sequences with P_A(s) + P_B(s) = -2 for s != 0.
Then H(2l + 2) exists (bordered two-circulant core; Fletcher-Gysin-Seberry).  For l = 2n - 1 this
gives H(4n).  Row sums are forced to be +-1 (r_A^2 + r_B^2 = 2), so WLOG both are +1.

Symmetry: H <= Z_l^*; sequences constant on H-orbits; condition on classes of <H, -1>.
Search: enumerate the whole pool (row sum +1), hash class vectors, join A with B. O(pool)."""
import sys, time
from math import gcd
import numpy as np

def subgroups(l, max_gens=2):
    units = [u for u in range(1, l) if gcd(u, l) == 1]
    def closure(gens):
        S, todo = {1}, [1]
        while todo:
            x = todo.pop()
            for g in gens:
                y = x * g % l
                if y not in S: S.add(y); todo.append(y)
        return frozenset(S)
    cyc = {}
    for u in units: cyc.setdefault(closure([u]), u)
    subs = {H: (g,) for H, g in cyc.items()}
    items = list(cyc.items())
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            subs.setdefault(closure([items[i][1], items[j][1]]), (items[i][1], items[j][1]))
    return subs

def structure(l, H):
    orb = -np.ones(l, dtype=np.int64); k = 0
    for t in range(l):
        if orb[t] < 0:
            for h in H: orb[t * h % l] = k
            k += 1
    L = set(H) | {(-h) % l for h in H}
    cls = -np.ones(l, dtype=np.int64); reps = []
    for s in range(1, l):
        if cls[s] < 0:
            for h in L: cls[s * h % l] = len(reps)
            reps.append(s)
    O, C = k, len(reps)
    S = np.zeros((C, O, O), dtype=np.int64)
    t = np.arange(l)
    for c, s in enumerate(reps):
        np.add.at(S[c], (orb, orb[(t + s) % l]), 1)
    sizes = np.bincount(orb)
    return orb, sizes, S, reps

def pool(sizes, r, chunk=1 << 20):
    O = len(sizes)
    out = []
    for start in range(0, 1 << O, chunk):
        idx = np.arange(start, min(start + chunk, 1 << O), dtype=np.int64)
        Y = 1 - 2 * ((idx[:, None] >> np.arange(O)) & 1)
        keep = (Y @ sizes) == r
        out.append(Y[keep].astype(np.int8))
    return np.concatenate(out) if out else np.zeros((0, O), np.int8)

def class_vectors(Y, S):
    Yf = Y.astype(np.float64)
    return np.rint(np.einsum("bo,cop,bp->bc", Yf, S.astype(np.float64), Yf, optimize=True)).astype(np.int64)

def find_pairs(l, H, maxpool=4e7, limit=5):
    orb, sizes, S, reps = structure(l, H)
    O = len(sizes)
    if 2 ** O > maxpool * 50: return None, O, 0
    Y = pool(sizes, 1)
    if len(Y) == 0 or len(Y) > maxpool: return None, O, len(Y)
    sols = []
    R = np.random.default_rng(1).integers(1, 2**62, size=S.shape[0], dtype=np.int64).astype(np.uint64) | np.uint64(1)
    target = np.full(S.shape[0], -2, dtype=np.int64)
    for s in range(0, len(Y), 1 << 18):
        V = class_vectors(Y[s:s + (1 << 18)], S)
        h = (V.astype(np.uint64) * R).sum(1, dtype=np.uint64)
        if s == 0: HV, VV = [h], [V]
        else: HV.append(h); VV.append(V)
    h = np.concatenate(HV); V = np.concatenate(VV)
    ht = (target.astype(np.uint64) * R).sum(dtype=np.uint64)
    order = np.argsort(h); hs = h[order]
    need = ht - h
    lo = np.searchsorted(hs, need, 'left'); hi = np.searchsorted(hs, need, 'right')
    for a in np.flatnonzero(hi > lo)[:1000]:
        for b in order[lo[a]:hi[a]]:
            if np.array_equal(V[a] + V[b], target):
                sols.append((Y[a], Y[b]))
                if len(sols) >= limit: return sols, O, len(Y)
    return sols, O, len(Y)

def to_sequences(l, H, sol):
    orb, sizes, S, reps = structure(l, H)
    return [np.array(y, dtype=np.int64)[orb] for y in sol]

def verify_lp(a, b):
    l = len(a)
    for s in range(1, l):
        if int(a @ np.roll(a, -s) + b @ np.roll(b, -s)) != -2: return False
    return True

if __name__ == "__main__":
    for l in map(int, sys.argv[1:]):
        t = time.time(); hits = 0
        for H, gens in sorted(subgroups(l).items(), key=lambda kv: -len(kv[0])):
            orb_count = None
            sols, O, N = find_pairs(l, H)
            if sols is None: continue
            if sols:
                a, b = to_sequences(l, H, sols[0])
                print(f"l={l}: |H|={len(H)} gens={gens} orbits={O} pool={N}: {len(sols)} pair(s); verified={verify_lp(a, b)}", flush=True)
                hits += 1
        print(f"l={l}: structures with a Legendre pair: {hits}  ({time.time()-t:.1f}s)", flush=True)

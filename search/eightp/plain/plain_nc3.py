"""plain_nc3.py -- plain multiplier structures for Spence over odd |G| = 2p-1 with G = Z_3^2 x Z_m (m coprime to 3) or G = Z_9 x Z_m, automorphism group
Aut(G) = GL(2,3) x Z_m^* (or Z_9^* x Z_m^*): random 2-generator subgroups, deduplicated by orbit partition, then the exact search of plain_spence3.
usage: plain_nc3.py KIND m nsamples [Omax] [maxL]   KIND = z3z3 | z9   (a random sample of subgroup shapes, not all subgroups: HEURISTIC coverage; each shape is solved exactly)"""
import sys, math, itertools, random, time
import numpy as np
import os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import plain_spence3 as ps

def gl23():
    els = []
    for a, b, c, d in itertools.product(range(3), repeat=4):
        if (a * d - b * c) % 3: els.append((a, b, c, d))
    return els

class Grp:
    def __init__(self, kind, m):
        self.kind = kind; self.m = m
        if kind == 'z3z3':
            self.v = 9 * m
            self.coords = np.array([(a1, a2, k) for a1 in range(3) for a2 in range(3) for k in range(m)], dtype=np.int64)
            self.mods = np.array([3, 3, m])
            self.A = gl23()
        else:
            self.v = 9 * m
            self.coords = np.array([(a, k) for a in range(9) for k in range(m)], dtype=np.int64)
            self.mods = np.array([9, m])
            self.A = [(u,) for u in range(1, 9) if math.gcd(u, 9) == 1]
        self.Um = [u for u in range(1, m) if math.gcd(u, m) == 1]
        self.strides = np.array([int(np.prod(self.mods[i + 1:])) for i in range(len(self.mods))], dtype=np.int64)
    def index(self, C): return (np.asarray(C) % self.mods) @ self.strides
    def act(self, el):
        """permutation of G induced by el = (Mtuple, u): (a, k) -> (M a, u k)"""
        M, u = el
        C = self.coords.copy()
        if self.kind == 'z3z3':
            a = C[:, :2]; Ma = np.stack([(M[0] * a[:, 0] + M[1] * a[:, 1]) % 3, (M[2] * a[:, 0] + M[3] * a[:, 1]) % 3], 1)
            C = np.concatenate([Ma, ((u * C[:, 2]) % self.m)[:, None]], 1)
        else:
            C = np.stack([(M[0] * C[:, 0]) % 9, (u * C[:, 1]) % self.m], 1)
        return self.index(C)
    def neg(self): return self.index(-self.coords)
    def shift(self, g): return self.index(self.coords + self.coords[g])

def orbits_of(perms, v):
    parent = list(range(v))
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for p in perms:
        for i in range(v):
            a, b = find(i), find(int(p[i]))
            if a != b: parent[a] = b
    lab = -np.ones(v, dtype=np.int64); reps = []
    roots = {}
    for i in range(v):
        r = find(i)
        if r not in roots: roots[r] = len(reps); reps.append(i)
        lab[i] = roots[r]
    return lab, reps

def run(kind, m, nsamp, Omax=28, maxL=2e5, seed=1):
    G = Grp(kind, m); v = G.v
    rng = random.Random(seed)
    zero = int(G.index(np.array([[0] * len(G.mods)]))[0]); assert zero == 0
    negp = G.neg()
    perm_cache = {}
    def P(el):
        if el not in perm_cache: perm_cache[el] = G.act(el)
        return perm_cache[el]
    seen = {}
    t0 = time.time()
    for s in range(nsamp):
        gens = []
        for _ in range(rng.choice([1, 2, 2, 3])):
            Mt = rng.choice(G.A); u = rng.choice(G.Um)
            gens.append((Mt, u))
        lab, reps = orbits_of([P(g) for g in gens], v)
        O = len(reps)
        if O > Omax or O < 3: continue
        key = hash(tuple(lab.tolist()))
        if key in seen: continue
        seen[key] = (lab, reps, gens)
    print(f'{kind} m={m} v={v}: {len(seen)} distinct partitions with 3..{Omax} orbits from {nsamp} samples ({time.time()-t0:.0f}s)', flush=True)
    res = []
    for key, (lab, reps, gens) in sorted(seen.items(), key=lambda kv: len(kv[1][1])):
        O = len(reps); size = np.bincount(lab, minlength=O)
        # row-sum feasibility
        tot = int(size.sum()); off = tot
        dp = np.zeros(2 * tot + 1); dp[off] = 1
        for sz in size:
            sz = int(sz); new = np.zeros_like(dp); new[sz:] += dp[:-sz]; new[:-sz] += dp[sz:]; dp = new
        L = dp[off + 1]
        if L == 0: continue
        negl = lab[negp[np.array(reps)]]
        sn = set(); cls = []
        for a in range(1, O):
            if a in sn: continue
            sn.add(a); sn.add(int(negl[a])); cls.append(a)
        rec = dict(kind=kind, m=m, O=O, classes=len(cls), L=float(L), gens=[(g[0], g[1]) for g in gens])
        if L > maxL: rec['status'] = 'bigL'; print(rec, flush=True); res.append(rec); continue
        Ms = []
        for a in cls:
            sh = G.shift(reps[a]); Mx = np.zeros((O, O), dtype=np.int64); np.add.at(Mx, (lab, lab[sh]), 1); Ms.append(Mx)
        X = ps.enum_list(O, size)
        Xi = X.astype(np.int64)
        V = np.stack([np.einsum('ni,ij,nj->n', Xi, Mx, Xi) for Mx in Ms], 1)
        Vu, inv = np.unique(V, axis=0, return_inverse=True); cnt = np.bincount(inv.ravel())
        tgt = np.full(len(cls), -4, dtype=np.int64)
        tt = time.time()
        tot_, found = ps.solve(Vu, cnt, tgt)
        rec.update(Lu=len(Vu), count_rho_plus=int(tot_), status='EMPTY' if tot_ == 0 else 'SOLUTION', secs=round(time.time() - tt, 1))
        if tot_:
            rec['found'] = [[X[np.flatnonzero(inv.ravel() == q)[0]].tolist() for q in quad] for quad in found[:1]]
            rec['orbit_labels'] = lab.tolist()
            print('FOUND', rec, flush=True)
        print({k: val for k, val in rec.items() if k not in ('found', 'orbit_labels')}, flush=True)
        res.append(rec)
    return res

if __name__ == '__main__':
    kind = sys.argv[1]; m = int(sys.argv[2]); n = int(sys.argv[3])
    Omax = int(sys.argv[4]) if len(sys.argv) > 4 else 24
    maxL = float(sys.argv[5]) if len(sys.argv) > 5 else 2e5
    run(kind, m, n, Omax, maxL)

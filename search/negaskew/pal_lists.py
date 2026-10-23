"""Paley-pattern exhaustive search for the quaternary skew problem of README.md on Z_v = Z_K x F_q (q = 3 mod 4 prime, K cyclic
of odd order, gcd(K, q) = 1), H = {1} x QR(q) x H_K, with at most 5 K-orbits (kappa <= 5).  For each defect-pattern class
(EXACT identity: sum over blocks X and k in K of [(q-1)|e|^2 + 2 Re(x0 conj e)] = 4|K|; A has a forced half defect at
k = 0) the explicit block lists are generated and joined exactly with qjoin.c.  Classes are reduced by the units of Z_K,
by x1 <-> x2 (nu) and by U -> conj U(-.) (sigma); a class whose lists have no row-sum combination is EXACTLY empty.
usage: python3 pal_lists.py v K q "[gens]" dmax [tag [part nparts]]
  gens  generators of H in Z_v (qshape.cyc_gens(K, q, H_K) gives them),  dmax  largest number of defects per block (5)
  environment: PAL_DRY=1 lists the classes, list sizes and the row-sum-aware cost (0 = no row-sum combination = EXACTLY
               empty) without joining;  PAL_CLASSES=1,5 only those classes;  PAL_NOSYM=1 no symmetry halving (full count);
               PAL_MAXPRINT=n solutions printed per class (default 3)
Needs ./qjoin (gcc -O3 -march=haswell -o qjoin qjoin.c).  Results are appended to runs/pal_<tag>.done.jsonl (a rerun
skips the finished classes); the last line is the total count over all classes.
Example (validation, 1m49s): PAL_NOSYM=1 python3 pal_lists.py 21 3 7 "[16]" 3 t21   ->   TOTAL ... 79164864 (= qskew over Z_21; without PAL_NOSYM the halved classes give 33882408)"""
import sys, json, subprocess, itertools, os, time
import numpy as np
from math import gcd
from collections import Counter, defaultdict
from qdrv import shape
from enum_shapes import gen_sub
from idfeas import korbits
from alloc import list_patterns

def units():
    return [(1, 0), (0, 1), (-1, 0), (0, -1)]
U = units()

def triple_info(q, e0, e1, e2):
    x0, x1, x2 = U[e0], U[e1], U[e2]
    d = (e2 - e1) & 3
    kind = 'P' if d == 2 else ('F' if d == 0 else 'M')
    es = (x1[0] + x2[0], x1[1] + x2[1])          # 2e
    n2 = es[0] ** 2 + es[1] ** 2                  # 4|e|^2
    re = x0[0] * es[0] + x0[1] * es[1]            # 2 Re(x0 conj e)... (x0 conj(2e)) real part = 2*Re(x0 conj e)*... ; Re(x0 conj es) = 2 Re(x0 conj e)
    W4 = (q - 1) * n2 + 4 * re                    # 4 * [ (q-1)|e|^2 + 2 Re(x0 conj e) ]
    return kind, W4

class PalShape:
    def __init__(self, v, K, q, gens):
        self.v, self.K, self.q = v, K, q
        assert gcd(K, q) == 1 and K * q == v
        self.d = shape(v, gens)
        d = self.d
        H = gen_sub(v, gens)
        HK = sorted(set(h % K for h in H))
        self.orbsK, self.idxK = korbits(K, HK)
        self.kap = len(self.orbsK)
        qr = set((x * x) % q for x in range(1, q))
        def typ(f): return 0 if f == 0 else (1 if f in qr else 2)
        self.orbid = {}
        for g in range(v):
            key = (self.idxK[g % K], typ(g % q))
            o = d['orb'][g]
            if key in self.orbid: assert self.orbid[key] == o, "orbit structure mismatch"
            else: self.orbid[key] = o
        assert len(self.orbid) == d['n'] == 3 * self.kap, (len(self.orbid), d['n'], self.kap)
        self.negK = [self.idxK[(-o[0]) % K] for o in self.orbsK]
        for (k, t), o in self.orbid.items():
            kp = self.negK[k]; tp = t if t == 0 else 3 - t
            assert d['neg'][o] == self.orbid[(kp, tp)]
        self.tnames = {}
        # option sets by type name
        self.opts = defaultdict(list)       # (kind,W4) -> list of (e0,e1,e2)
        for e0 in range(4):
            for e1 in range(4):
                for e2 in range(4):
                    kind, W4 = triple_info(q, e0, e1, e2)
                    self.opts[(kind, W4)].append((e0, e1, e2))
    def type_key(self, name):
        q = self.q
        if name == 'F+': return ('F', 4 * (q + 1))
        if name == 'F-': return ('F', 4 * (q - 3))
        if name == 'Fi': return ('F', 4 * (q - 1))
        if name == 'M+': return ('M', 2 * (q + 1))
        if name == 'M-': return ('M', 2 * (q - 3))
        raise ValueError(name)

    def block_list(self, cfg, kind):
        """kind 'B' (free block, e at orbit (0,0) fixed to 0) or ('A', w0name).  cfg = tuple of (korbit, typename) defects.
        returns np.uint64 array of packed exponent vectors"""
        kap = self.kap; n = self.d['n']; opts = self.opts
        dmap = dict(cfg)
        if kind == 'B':
            per = []
            for i in range(kap):
                if i in dmap: o = opts[self.type_key(dmap[i])]
                else: o = opts[('P', 0)]
                if i == 0: o = [t for t in o if t[0] == 0]
                per.append(o)
            return self._pack_free(per)
        else:
            w0name = kind[1]
            # kappa = 0 : forced mixed ; e0 in {0,3}; e2 = 1 - e1
            o0 = []
            tk = self.type_key(w0name)
            for e0 in (0, 3):
                for e1 in range(4):
                    e2 = (1 - e1) & 3
                    if triple_info(self.q, e0, e1, e2) == tk: o0.append((e0, e1, e2))
            pairs = [i for i in range(1, kap) if i < self.negK[i]]
            assert all(self.negK[i] != i for i in range(1, kap))
            per = []
            for i in pairs:
                if i in dmap: o = opts[self.type_key(dmap[i])]
                else: o = opts[('P', 0)]
                per.append(o)
            return self._pack_A(o0, pairs, per)

    def _pack_free(self, per):
        n = self.d['n']; sizes = [len(o) for o in per]
        tot = int(np.prod(sizes, dtype=object))
        if tot == 0: return np.zeros(0, dtype=np.uint64)
        idx = np.arange(tot, dtype=np.int64)
        out = np.zeros(tot, dtype=np.uint64)
        for i, o in enumerate(per):
            m = len(o); sel = idx % m; idx = idx // m
            arr = np.array(o, dtype=np.uint64)   # m x 3
            for t in range(3):
                oid = self.orbid[(i, t)]
                out |= (arr[sel, t] << np.uint64(2 * oid))
        return out

    def _pack_A(self, o0, pairs, per):
        n = self.d['n']
        sizes = [len(o0)] + [len(o) for o in per]
        tot = int(np.prod(sizes, dtype=object))
        if tot == 0: return np.zeros(0, dtype=np.uint64)
        idx = np.arange(tot, dtype=np.int64)
        out = np.zeros(tot, dtype=np.uint64)
        # kappa 0
        m = len(o0); sel = idx % m; idx = idx // m; arr = np.array(o0, dtype=np.uint64)
        for t in range(3):
            out |= (arr[sel, t] << np.uint64(2 * self.orbid[(0, t)]))
        for j, i in enumerate(pairs):
            o = per[j]; m = len(o); sel = idx % m; idx = idx // m; arr = np.array(o, dtype=np.uint64)
            ip = self.negK[i]
            for t in range(3):
                out |= (arr[sel, t] << np.uint64(2 * self.orbid[(i, t)]))
            # partner: neg orbits: e'[(ip,0)] = 1 - e0 ; e'[(ip,2)] = 1 - e1 ; e'[(ip,1)] = 1 - e2
            for (src, dst) in ((0, 0), (1, 2), (2, 1)):
                val = (np.uint64(1) - arr[sel, src]) & np.uint64(3)
                val = (np.uint64(5) - arr[sel, src]) & np.uint64(3)   # 1 - x mod 4 computed as (5 - x) mod 4 to stay unsigned
                out |= (val << np.uint64(2 * self.orbid[(ip, dst)]))
        return out

def unit_canon(ps, key):
    import math
    (cA, w0), cB, cC, cD = key
    K = ps.K; best = None
    for u in range(1, K):
        if math.gcd(u, K) != 1: continue
        pi = [ps.idxK[(u * ps.orbsK[k][0]) % K] for k in range(ps.kap)]
        def mp(cfg, paired=False):
            out = []
            for (k, t) in cfg:
                k2 = pi[k]
                if paired: k2 = min(k2, ps.negK[k2])
                out.append((k2, t))
            return tuple(sorted(out))
        cand = ((mp(cA, True), w0),) + tuple(sorted([mp(cB), mp(cC), mp(cD)]))
        if best is None or str(cand) < str(best): best = cand
    return best

def swap_ex(ps, arr):
    out = arr.copy()
    for i in range(ps.kap):
        o1 = ps.orbid[(i, 1)]; o2 = ps.orbid[(i, 2)]
        m1 = np.uint64(3) << np.uint64(2 * o1); m2 = np.uint64(3) << np.uint64(2 * o2)
        v1 = (arr >> np.uint64(2 * o1)) & np.uint64(3); v2 = (arr >> np.uint64(2 * o2)) & np.uint64(3)
        out = (out & ~(m1 | m2)) | (v2 << np.uint64(2 * o1)) | (v1 << np.uint64(2 * o2))
    return out

def sigma_ex(ps, arr):
    # U'(g) = conj U(-g): orbit (i,t) -> (negK[i], t' = t or 3-t), exponent e -> -e
    out = np.zeros_like(arr)
    for i in range(ps.kap):
        for t in range(3):
            src = ps.orbid[(i, t)]; dst = ps.orbid[(ps.negK[i], t if t == 0 else 3 - t)]
            val = (np.uint64(4) - ((arr >> np.uint64(2 * src)) & np.uint64(3))) & np.uint64(3)
            out |= (val << np.uint64(2 * dst))
    return out

def canon_classes(pats):
    cn = Counter(); rep = {}
    for (a, b, c, d) in pats:
        key = (a,) + tuple(sorted([b, c, d]))
        cn[key] += 1
    return cn

def run_class(ps, key, mult, tag, workdir='runs', maxprint=3, dry=False):
    (cA, w0), cB, cC, cD = key
    q = ps.q
    w0name = 'M+' if w0 == (q + 1) // 2 else 'M-'
    L = [ps.block_list(cA, ('A', w0name)), ps.block_list(cB, 'B'), ps.block_list(cC, 'B'), ps.block_list(cD, 'B')]
    # symmetries (all preserve PAF equations; see STATUS): sigma_X : U -> conj U(-.) acts on each free block independently (keeps the pattern iff all its
    # defects sit at k = 0); nu : x1 <-> x2 on every orbit of all blocks.  Use sigma on every eligible free block; else nu on the largest list.
    halved = []
    elig = [all(k == 0 for (k, t) in cfg) for cfg in (cB, cC, cD)]
    if os.environ.get('PAL_NOSYM'): elig = [False, False, False]; nu_ok = False
    else: nu_ok = True
    if any(elig):
        for j in range(3):
            if elig[j]:
                sg = sigma_ex(ps, L[j + 1]); L[j + 1] = L[j + 1][L[j + 1] <= sg]; halved.append(('sigma', j + 1))
    elif nu_ok:
        big = max((1, 2, 3), key=lambda i: len(L[i]))
        sw = swap_ex(ps, L[big]); L[big] = L[big][L[big] <= sw]; halved.append(('nu', big))
    sizes = [len(x) for x in L]
    # DC-aware cost of each ordered pairing (table = first pair, probes = second pair)
    n_ = ps.d['n']; szv = np.array(ps.d['sz']); v_ = ps.v; tot = 4 * v_
    UXa = np.array([1, 0, -1, 0]); UYa = np.array([0, 1, 0, -1])
    def Rh(arr):
        rr = np.zeros(len(arr), dtype=np.int64); ri = np.zeros(len(arr), dtype=np.int64)
        for a in range(n_):
            e = ((arr >> np.uint64(2 * a)) & np.uint64(3)).astype(np.int64)
            rr += szv[a] * UXa[e]; ri += szv[a] * UYa[e]
        return Counter(zip(rr.tolist(), ri.tolist()))
    Hh = [Rh(x) for x in L]
    def pairstat(i, j):
        el = Counter(); bk = Counter()
        for (a, b), c in Hh[i].items():
            na = a * a + b * b
            for (x, y), dd in Hh[j].items():
                sm = na + x * x + y * y
                if sm <= tot: el[sm] += c * dd; bk[sm] += 1
        return el, bk
    res = {}
    for (p, r) in (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))):
        e1, b1 = pairstat(*p); e2, b2 = pairstat(*r)
        res[(p, r)] = (sum(e1[s_] for s_ in e1 if (tot - s_) in e2), sum(b1[s_] * e2.get(tot - s_, 0) for s_ in b1))
        res[(r, p)] = (sum(e2[s_] for s_ in e2 if (tot - s_) in e1), sum(b2[s_] * e1.get(tot - s_, 0) for s_ in b2))
    bestk = min(res, key=lambda k: res[k][0] * 3 + res[k][1])
    cost = res[bestk][0] * 3 + res[bestk][1]
    p, r = bestk
    if cost == 0:
        return sizes, (p, r), None, halved, cost
    fn = f"{workdir}/pl_{tag}"
    for i in range(4): L[i].tofile(f"{fn}.l{i}.bin")
    if dry:
        return dict(sizes=sizes, pairing=(p, r), cost=cost)
    return sizes, (p, r), fn, halved, cost

if __name__ == "__main__":
    v = int(sys.argv[1]); K = int(sys.argv[2]); q = int(sys.argv[3]); gens = json.loads(sys.argv[4]); dmax = int(sys.argv[5])
    tag = sys.argv[6] if len(sys.argv) > 6 else f"{v}"
    part = int(sys.argv[7]) if len(sys.argv) > 7 else 0; nparts = int(sys.argv[8]) if len(sys.argv) > 8 else 1
    os.makedirs('shapes', exist_ok=True); os.makedirs('runs', exist_ok=True)
    donef = f"runs/pal_{tag}.done.jsonl"; done = set()
    if os.path.exists(donef):
        for l in open(donef): done.add(json.loads(l)['cls'])
    ps = PalShape(v, K, q, gens)
    shp = f"shapes/pal_{tag}.txt"; shape(v, gens, shp)
    pats, szK = list_patterns(v, K, q, gens, dmax)
    cn0 = canon_classes(pats)
    cn = Counter()
    for k0, c0 in cn0.items(): cn[unit_canon(ps, k0)] += c0
    print('ordered patterns', len(pats), 'classes', len(cn), flush=True)
    keys = sorted(cn.keys(), key=lambda k: str(k))
    total = 0; results = []
    for ci, key in enumerate(keys):
        if ci % nparts != part or ci in done: continue
        if os.environ.get('PAL_CLASSES') and str(ci) not in os.environ['PAL_CLASSES'].split(','): continue
        t0 = time.time()
        r = run_class(ps, key, cn[key], f"{tag}_{ci}")
        sizes, (p, rr), fn, halved, cost_ = r
        if os.environ.get('PAL_DRY'):
            print(ci, cn[key], str(key), sizes, 'cost3b+p %.3g' % cost_, halved, flush=True)
            if fn is not None:
                for i in range(4):
                    try: os.remove(f"{fn}.l{i}.bin")
                    except OSError: pass
            continue
        if fn is None:
            rec = dict(cls=ci, key=str(key), mult=cn[key], sizes=sizes, count=0, note='no DC-valid combination', sec=round(time.time() - t0, 1)); results.append(rec); print(json.dumps(rec), flush=True); open(donef, 'a').write(json.dumps(rec) + '\n'); continue
        cmd = ['./qjoin', shp, f"{fn}.l0.bin", f"{fn}.l1.bin", f"{fn}.l2.bin", f"{fn}.l3.bin", str(p[0]), str(p[1]), str(rr[0]), str(rr[1]), os.environ.get('PAL_MAXPRINT','3'), '15000000']
        res = subprocess.run(cmd, capture_output=True, text=True)
        out = res.stdout
        cnt = None
        for l in out.split('\n'):
            if l.startswith('COUNT'): cnt = int(l.split()[1])
        sols = [l for l in out.split('\n') if l.startswith('SOL')]
        for i in range(4):
            try: os.remove(f"{fn}.l{i}.bin")
            except OSError: pass
        total += (cnt or 0) * cn[key]
        rec = dict(cls=ci, key=str(key), mult=cn[key], sizes=sizes, pairing=[list(p), list(rr)], count=cnt, sec=round(time.time() - t0, 1), sols=sols[:3] if not os.environ.get('PAL_MAXPRINT') else sols, halved=halved)
        results.append(rec); print(json.dumps(rec), flush=True); open(donef, 'a').write(json.dumps(rec) + '\n')
    print('TOTAL ordered-pattern count (normalised, A skew list, B C D with e0 = 0):', total)

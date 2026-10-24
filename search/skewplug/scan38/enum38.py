"""Enumerate EVERY real K-symmetric census shape (skew Paley plug-in, G = K x F_q, q = 3 mod 4 prime power, H = H_K x squares).
v in search/data/skew_census_undecided.txt and not settled in the repository (skew/N with N = 4 mod 8, or hadamard/3732 = v 933, 919), v = |K| q, |K| odd >= 3.
K: every abelian group of order |K| (cyclic; non-cyclic by Sylow types).  H_K: odd-order subgroups of Aut(K):
  cyclic K                : ALL odd-order subgroups of the unit group of Z_m
  Z_p^2 blocks            : ALL odd-order subgroups of GL(2,p) up to conjugacy (diagonal subgroups of the split torus, cyclic subgroups of the Singer cycle, U x| D with U unipotent)
  Z_p^k (k >= 3) blocks   : diagonal-torus and Singer cyclic subgroups  (partial coverage; flagged)
  other non-cyclic parts  : diagonal units per component (partial; flagged)
Output shapes38.json: list of dicts (v, q, kmods, hdesc, order, nu, ...).   usage (in search/skewplug): python3 scan38/enum38.py [NUMAX=25] [OUT=shapes38.json]   (the list depends on skew/: values settled since the scan, such as 385, are not enumerated)
"""
import sys, os, json, itertools, math, time, warnings
from math import gcd
warnings.filterwarnings('ignore')
import numpy as np
from sympy import factorint
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..'))      # scan38/ and the modules of search/skewplug
DATA = os.path.join(HERE, 'data')
REPO = os.environ.get('REPO', os.path.join(HERE, '..', '..', '..'))      # the repository root (search/skewplug/scan38/ is three levels down)

# ---------------------------------------------------------------- census targets
def targets():
    und = [int(l) for l in open(REPO + '/search/data/skew_census_undecided.txt') if l.strip() and not l.startswith('#')]
    sk = [int(d) for d in os.listdir(REPO + '/skew') if d.isdigit()]
    settled = {N // 4 for N in sk if N % 8 == 4} | {933, 919}      # 933: hadamard/3732, 919: settled by the repository (r37 list)
    return [v for v in und if v not in settled]

def qlist(v):
    out = []
    for p, e in factorint(v).items():
        if p % 4 == 3:
            for k in range(1, e + 1):
                q = p ** k
                if q % 4 == 3: out.append(q)
    return sorted(out)

def identity_ok(m, q):
    """coarse: exists N_w >= 1, N_- with (q+1) N_w - 4 N_- = 4m, 0 <= N_- <= N_w"""
    if q == 3: return True
    return any(4 * m >= Nw * (q - 3) and 4 * m <= Nw * (q + 1) and ((q + 1) * Nw - 4 * m) % 4 == 0 for Nw in range(1, 4 * m // (q - 3) + 2))

# ---------------------------------------------------------------- group elements (tuples of per-block actions)
# blocks: ('cyc', n) -> int unit; ('mat', p, k) -> tuple of k rows (each a tuple), acting on column vectors
def bl_mul(b, x, y):
    if b[0] == 'cyc': return (x * y) % b[1]
    p, k = b[1], b[2]
    return tuple(tuple(sum(x[i][t] * y[t][j] for t in range(k)) % p for j in range(k)) for i in range(k))
def bl_id(b):
    if b[0] == 'cyc': return 1
    k = b[2]; return tuple(tuple(1 if i == j else 0 for j in range(k)) for i in range(k))
def el_mul(blocks, x, y): return tuple(bl_mul(b, a, c) for b, a, c in zip(blocks, x, y))
def el_id(blocks): return tuple(bl_id(b) for b in blocks)
def closure(blocks, gens, cap=5000):
    e = el_id(blocks); G = {e}; fr = [e]
    while fr:
        x = fr.pop()
        for g in gens:
            y = el_mul(blocks, x, g)
            if y not in G:
                G.add(y); fr.append(y)
                if len(G) > cap: return None
    return frozenset(G)
def el_order(blocks, g, cap=100000):
    e = el_id(blocks); x = g; o = 1
    while x != e:
        x = el_mul(blocks, x, g); o += 1
        if o > cap: return None
    return o

def flat_mods(blocks):
    out = []
    for b in blocks:
        out += [b[1]] if b[0] == 'cyc' else [b[1]] * b[2]
    return out

def perm_of(blocks, el, el_coords, w):
    """permutation of K indices induced by group element el"""
    cols = []; c = 0
    for b, a in zip(blocks, el):
        if b[0] == 'cyc':
            cols.append((el_coords[:, c] * a) % b[1]); c += 1
        else:
            p, k = b[1], b[2]; M = np.array(a, dtype=np.int64)
            X = el_coords[:, c:c + k]
            Y = (X @ M.T) % p
            for j in range(k): cols.append(Y[:, j])
            c += k
    new = np.stack(cols, 1)
    return (new * w[None, :]).sum(1)

# ---------------------------------------------------------------- per-block odd subgroup lists
def odd_units(n):
    out = []
    for u in range(1, n):
        if gcd(u, n) != 1: continue
        x, o = u, 1
        while x != 1: x = x * u % n; o += 1
        if o % 2 == 1: out.append(u)
    return out
def subgroups_abelian(blocks, elems, cap_groups=3000):
    """all subgroups of the abelian group formed by elems (closed under multiplication, containing identity): cyclic subgroups, then closure under joins"""
    e = el_id(blocks); start = frozenset([e])
    cyc = {}
    for g in elems:
        if g == e: continue
        C = closure(blocks, [g])
        if C is not None: cyc[C] = g
    cl = list(cyc.items())
    subs = {start}; fr = [start]
    for C, g in cl: 
        if C not in subs: subs.add(C); fr.append(C)
    while fr:
        S = fr.pop()
        for C, g in cl:
            if C <= S: continue
            T = closure(blocks, list(S) + [g]) if len(S) < 200 else frozenset(el_mul(blocks, x, y) for x in S for y in C)
            if T is not None and T not in subs:
                subs.add(T); fr.append(T)
                if len(subs) > cap_groups: return subs
    return subs

def cyc_block_subgroups(n):
    b = [('cyc', n)]
    el = [(u,) for u in odd_units(n)]
    # the odd elements of the unit group form a subgroup (odd part)
    return b, subgroups_abelian(b, el + [(1,)])

def mat2_subgroups(p):
    """all odd-order subgroups of GL(2,p) up to conjugacy: (a) subgroups of the split torus (odd part), (b) cyclic subgroups of the Singer cycle not in (a), (c) U x| D, D <= T"""
    blk = [('mat', p, 2)]
    ou = odd_units(p)
    T = [(((a, 0), (0, b)),) for a in ou + [1] for b in ou + [1]]
    T = sorted(set(T))
    subsT = subgroups_abelian(blk, T)
    out = {}
    for S in subsT: out[S] = 'T'
    # (c) U x| D
    Umat = (((1, 1), (0, 1)),)
    for S in subsT:
        H = closure(blk, list(S) + [Umat])
        if H is not None and H not in out: out[H] = 'UD'
    # (b) Singer: cyclic subgroups of F_{p^2}^* of odd order d not dividing p-1
    c = next(c for c in range(2, p) if pow(c, (p - 1) // 2, p) == p - 1) if p > 2 else 1
    def mul(u, v): return ((u[0] * v[0] + c * u[1] * v[1]) % p, (u[0] * v[1] + u[1] * v[0]) % p)
    def power(u, n):
        r = (1, 0)
        while n:
            if n & 1: r = mul(r, u)
            u = mul(u, u); n >>= 1
        return r
    full = p * p - 1
    def prim():
        pf = [f for f in factorint(full)]
        for a in range(p):
            for b in range(1, p):
                g = (a, b)
                if all(power(g, full // f) != (1, 0) for f in pf): return g
    g = prim()
    odd_part = full
    while odd_part % 2 == 0: odd_part //= 2
    for d in sorted(d for d in range(1, odd_part + 1) if odd_part % d == 0):
        if d == 1 or (p - 1) % d == 0: continue      # scalars are in T
        h = power(g, full // d)
        Mx = (((h[0], (c * h[1]) % p), (h[1], h[0])),)
        H = closure(blk, [Mx[0:1][0:1][0] if False else Mx])
        if H is not None: out.setdefault(H, 'S')
    return blk, out

def _matpow(M, e, p):
    k = M.shape[0]; R = np.eye(k, dtype=np.int64)
    while e:
        if e & 1: R = (R @ M) % p
        M = (M @ M) % p; e >>= 1
    return R

def matk_subgroups(p, k):
    """Z_p^k (k >= 3): odd subgroups of the diagonal torus (scalars and diagonal), plus cyclic Singer subgroups of F_{p^k}^* of odd order (companion matrix of a primitive polynomial); partial coverage"""
    blk = [('mat', p, k)]
    ou = odd_units(p) + [1]
    T = sorted(set((tuple(tuple(d[i] if i == j else 0 for j in range(k)) for i in range(k)),) for d in itertools.product(ou, repeat=k)))
    subs = subgroups_abelian(blk, T, cap_groups=400)
    out = {S: 'T' for S in subs}
    qk = p ** k; full = qk - 1
    pf = list(factorint(full))
    I = np.eye(k, dtype=np.int64); C = None
    for coeffs in itertools.product(range(p), repeat=k):
        if coeffs[0] == 0: continue
        M = np.zeros((k, k), dtype=np.int64)
        for j in range(k - 1): M[j + 1, j] = 1
        for i in range(k): M[i, k - 1] = (-coeffs[i]) % p
        if not (_matpow(M, full, p) == I).all(): continue
        if all(not (_matpow(M, full // f, p) == I).all() for f in pf): C = M; break
    assert C is not None
    odd_part = full
    while odd_part % 2 == 0: odd_part //= 2
    for d in sorted(d for d in range(2, odd_part + 1) if odd_part % d == 0):
        g = _matpow(C, full // d, p)
        Mx = tuple(tuple(int(x) for x in row) for row in g)
        H = closure(blk, [(Mx,)], cap=3000)
        if H is not None: out.setdefault(H, 'S')
    return blk, out

# ---------------------------------------------------------------- K types
def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0: yield (); return
    for k in range(min(n, maxpart), 0, -1):
        for p in partitions(n - k, k): yield (k,) + p

def k_types(m):
    """abelian groups of order m (odd).  Returns list of (label, blocks_spec) where blocks_spec = list of block descriptors:
       ('cyc', n)  : cyclic part Z_n with unit actions (all coprime cyclic p-parts merged),
       ('mat', p, k): Z_p^k elementary block with matrix actions (k >= 2),
       diag components of non-elementary p-parts: ('cyc', p^a) separate blocks"""
    f = factorint(m); per = []
    for p, e in f.items():
        opts = []
        for part in partitions(e):
            opts.append((p, part))
        per.append(opts)
    res = []
    for combo in itertools.product(*per):
        cyc_n = 1; blocks = []; extra = []; label = []
        for (p, part) in combo:
            label.append('%d^(%s)' % (p, ','.join(map(str, part))))
            if len(part) == 1: cyc_n *= p ** part[0]
            elif all(a == 1 for a in part): blocks.append(('mat', p, len(part)))
            else:
                # largest component merged into the cyclic block? no: components of non-elementary parts as separate cyc blocks (diag actions only)
                for a in part: extra.append(('cyc', p ** a))
        bl = ([('cyc', cyc_n)] if cyc_n > 1 else []) + blocks + extra
        res.append((' x '.join(label), bl, any(len(part) > 1 and not all(a == 1 for a in part) for (p, part) in combo)))
    return res

# ---------------------------------------------------------------- shape invariants
def wl_signature(perm_list, N, w, negperm=None):
    """orbits of the group generated by perm_list on K (indices) and a Weisfeiler-Leman-type signature of the orbit structure under addition"""
    orb = -np.ones(N, dtype=np.int64); reps = []
    for k0 in range(N):
        if orb[k0] >= 0: continue
        o = len(reps); reps.append(k0); orb[k0] = o; st = [k0]
        while st:
            k = st.pop()
            for P in perm_list:
                k2 = int(P[k])
                if orb[k2] < 0: orb[k2] = o; st.append(k2)
    nu = len(reps); size = np.bincount(orb, minlength=nu)
    return orb, reps, nu, size

def add_table(mods, el, w):
    m = np.array(mods)
    def addk(k0):
        return (((el[k0][None, :] + el) % m) * w).sum(1)
    return addk

def signature(orb, reps, nu, size, addk, N):
    lab = [int(size[o]) for o in range(nu)]
    for rnd in range(3):
        new = []
        for o in range(nu):
            s = addk(reps[o])        # k0 + k' for all k'
            lo = np.array(lab)[orb]  # label of orbit of k'
            lo2 = np.array(lab)[orb[s]]
            from collections import Counter
            c = Counter(zip(lo.tolist(), lo2.tolist()))
            new.append(hash((lab[o], tuple(sorted(c.items())))) & 0xFFFFFFFFFFFF)
        lab = new
    return tuple(sorted(lab))

# ---------------------------------------------------------------- main enumeration
def enum_shapes(NUMAX=25, only=None, verbose=True):
    from sk import SShape
    from sk_pat import allocations, m_splits
    rows = []
    t0 = time.time()
    tg = targets()
    pairs = []
    for v in tg:
        for q in qlist(v):
            m = v // q
            if m < 3 or m % 2 == 0: continue
            if not identity_ok(m, q): continue
            pairs.append((v, q, m))
    if only: pairs = [p for p in pairs if p[0] in only]
    if verbose: print('census targets', len(tg), '(v, q) pairs passing the coarse identity', len(pairs), flush=True)
    for (v, q, m) in pairs:
        for (label, blocks, partial) in k_types(m):
            mods = flat_mods(blocks)
            N = int(np.prod(mods))
            if N != m: continue
            # per-block subgroup lists
            lists = []
            for b in blocks:
                if b[0] == 'cyc': lists.append(list(cyc_block_subgroups(b[1])[1]))
                elif b[2] == 2: lists.append(mat2_subgroups(b[1])[1].keys())
                else: lists.append(matk_subgroups(b[1], b[2])[1].keys())
            lists = [list(l) for l in lists]
            # subgroups of the product: direct products, and cyclic subdirect subgroups <(g1, ..., gk)>
            cand = set()
            ide = el_id(blocks)
            def prod_group(parts):
                gens = []
                for i, S in enumerate(parts):
                    for g in S:
                        t = list(ide); t[i] = g[0]; gens.append(tuple(t))
                return closure(blocks, gens, cap=20000)
            combos = list(itertools.product(*lists)) if np.prod([len(l) for l in lists], dtype=np.float64) < 5e4 else None
            if combos is None:
                if verbose: print('  skip (too many block-subgroup combos)', v, q, label, [len(l) for l in lists], flush=True)
                continue
            for parts in combos:
                H = prod_group(parts)
                if H is not None and len(H) % 2 == 1: cand.add(H)
            if len(blocks) > 1:
                # cyclic subdirect subgroups from element tuples: one generator per block taken from the block subgroups' generators
                gl = []
                for i, L in enumerate(lists):
                    gset = set()
                    for S in L:
                        for g in S: gset.add(g[0])
                    gl.append(sorted(gset, key=str)[:400])
                if np.prod([len(g) for g in gl], dtype=np.float64) < 2e5:
                    for tup in itertools.product(*gl):
                        H = closure(blocks, [tuple(tup)], cap=5000)
                        if H is not None and len(H) % 2 == 1: cand.add(H)
            el_coords = np.array(list(itertools.product(*[range(x) for x in mods])), dtype=np.int64).reshape(N, len(mods))
            wv = np.ones(len(mods), dtype=np.int64)
            for i in range(len(mods) - 2, -1, -1): wv[i] = wv[i + 1] * mods[i + 1]
            addk = add_table(mods, el_coords, wv)
            seen_sig = {}
            for H in sorted(cand, key=lambda S: (len(S), sorted(map(str, S)))):
                # gens: all elements is too many for orbit computation: use a small generating set found greedily
                gens = []; cur = frozenset([ide])
                for g in sorted(H, key=str):
                    if g in cur: continue
                    gens.append(g); cur = closure(blocks, gens)
                    if len(cur) == len(H): break
                perms = [perm_of(blocks, g, el_coords, wv) for g in gens]
                orb, reps, nu, size = wl_signature(perms, N, wv)
                if nu > NUMAX: continue
                sig = signature(orb, reps, nu, size, addk, N)
                key = (len(H), nu, sig)
                if key in seen_sig: continue
                seen_sig[key] = 1
                rows.append(dict(v=v, q=q, m=m, ktype=label, blocks=[list(b) for b in blocks], kmods=mods, order=len(H), nu=nu, partial=bool(partial),
                                 gens=[json.loads(json.dumps(g)) for g in gens], sizes=sorted(size.tolist()), sig=[str(x) for x in sig][:0]))
        if verbose: print('  done', v, q, m, 'rows so far', len(rows), round(time.time() - t0), 's', flush=True)
    return rows

if __name__ == '__main__':
    NUMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 25
    out = sys.argv[2] if len(sys.argv) > 2 else 'shapes38.json'
    only = set(int(x) for x in os.environ['ONLY'].split(',')) if os.environ.get('ONLY') else None
    rows = enum_shapes(NUMAX, only)
    json.dump(rows, open(out, 'w'))
    from collections import Counter
    print(len(rows), 'shapes (pre feasibility)', Counter(r['nu'] for r in rows))

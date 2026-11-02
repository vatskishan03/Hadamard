"""Shared helpers: orbit tables for (w, h), reading engine row files, expansion under the symmetry group, direct PAF checks."""
import os, struct, glob, itertools, subprocess, re

def is_prime(n):
    return n > 1 and all(n % d for d in range(2, int(n ** .5) + 1))

def primroot(p):
    for a in range(2, p):
        if all(pow(a, (p - 1) // q, p) != 1 for q in range(2, p) if (p - 1) % q == 0 and is_prime(q)):
            return a

class Params:
    def __init__(self, w, h):
        self.w, self.h = w, h
        self.m = (w - 1) // h
        self.h2 = self.m // 2
        self.g0 = primroot(w)
        self.dlog = {}
        x = 1
        for i in range(w - 1):
            self.dlog[x] = i
            x = x * self.g0 % w
        self.cls = {x: d % self.m for x, d in self.dlog.items()}

    def expand(self, xinf, orb):
        """orbit form -> list on Z_w"""
        return [xinf] + [orb[self.cls[x]] for x in range(1, self.w)]

def mask_to_orb(pos, neg, m):
    return [(1 if pos >> j & 1 else (-1 if neg >> j & 1 else 0)) for j in range(m)]

def s2(v):
    return "".join("+" if x > 0 else ("-" if x < 0 else "0") for x in v)

def read_rows(outdir, types):
    """engine rows -> list of (type_id, Xinf, Xorb(list), Yinf, Yorb)"""
    rows = []
    P = types["P"]
    for fn in glob.glob(os.path.join(outdir, "rows", "t*_k*.bin")):
        data = open(fn, "rb").read()
        for i in range(0, len(data), 20):
            xp, xn, yp, yn, t, _ = struct.unpack_from("<IIIIHH", data, i)
            ty = types[t]
            rows.append((t, ty["exi"], mask_to_orb(xp, xn, P.m), ty["eyi"], mask_to_orb(yp, yn, P.m)))
    return rows

def parse_types(w, h, swap_reduce=False, pgp="./src/pgp"):
    out = subprocess.run([pgp, "info", str(w), str(h)] + (["--swap-reduce"] if swap_reduce else []), capture_output=True, text=True).stdout
    types = {"P": Params(w, h)}
    for ln in out.splitlines():
        mm = re.match(r"type (\d+): sums \((\d+),(\d+)\)\s+X_inf=(-?\d+) sigma_X=(-?\d+) alpha_X=(\d+)\s+Y_inf=(-?\d+) sigma_Y=(-?\d+) alpha_Y=(\d+)", ln)
        if mm:
            t = int(mm.group(1))
            types[t] = dict(x=int(mm.group(2)), y=int(mm.group(3)), exi=int(mm.group(4)), sx=int(mm.group(5)), ax=int(mm.group(6)),
                            eyi=int(mm.group(7)), sy=int(mm.group(8)), ay=int(mm.group(9)))
    return types

def paf_cyc(seq):
    n = len(seq)
    return [sum(seq[i] * seq[(i + s) % n] for i in range(n)) for s in range(n)]


# ---------------------------------------------------------------- solutions: assembly, symmetry group, direct check
def seq_from_masks(pos, neg, m):
    return mask_to_orb(pos, neg, m)

def assemble(P, A, B, C, D, ainf, binf, cinf, dinf):
    """A..D orbit lists (length m), infinity values -> (a0,a1,b0,b1) as lists [inf, orbits...]; raises if not +-1"""
    def comb(X, xi, Y, yi, sgn):
        out = [xi + sgn * yi] + [X[j] + sgn * Y[j] for j in range(P.m)]
        if any(abs(v) != 1 for v in out):
            raise ValueError("supports are not complementary")
        return out
    return (comb(A, ainf, C, cinf, 1), comb(A, ainf, C, cinf, -1), comb(B, binf, D, dinf, 1), comb(B, binf, D, dinf, -1))

def rot(v, r, m):
    return [v[0]] + [v[1 + (j + r) % m] for j in range(m)]

def neg(v):
    return [-x for x in v]

def gen_group_images(sol, m):
    """all images of sol=(a0,a1,b0,b1) under: multiplier rotations, reversal of a alone / b alone, negations,
    translations of a / b by the element of order 2, swap a<->b, alternation (a1,b1 -> -a1,-b1)"""
    def key(s): return tuple(tuple(x) for x in s)
    gens = [
        lambda s: tuple(rot(x, 1, m) for x in s),
        lambda s: (rot(s[0], m // 2, m), rot(s[1], m // 2, m), s[2], s[3]),
        lambda s: (s[0], s[1], rot(s[2], m // 2, m), rot(s[3], m // 2, m)),
        lambda s: (neg(s[0]), neg(s[1]), s[2], s[3]),
        lambda s: (s[0], s[1], neg(s[2]), neg(s[3])),
        lambda s: (s[1], s[0], s[2], s[3]),
        lambda s: (s[0], s[1], s[3], s[2]),
        lambda s: (s[2], s[3], s[0], s[1]),
        lambda s: (s[0], neg(s[1]), s[2], neg(s[3])),
    ]
    seen = {key(sol): sol}
    todo = [sol]
    while todo:
        s = todo.pop()
        for g in gens:
            t = g(s); k = key(t)
            if k not in seen:
                seen[k] = t; todo.append(t)
    return seen

def sol_str(sol):
    return tuple(s2(x) for x in sol)

def pgp_check(P, sol):
    """exact check on the cyclic group Z_2w (CRT index y <-> (y mod 2, y mod w)) from the definition"""
    a0, a1, b0, b1 = sol
    w = P.w
    def full(l0, l1):
        out = []
        for y in range(2 * w):
            e, r = y % 2, y % w
            l = l1 if e else l0
            out.append(l[0] if r == 0 else l[1 + P.cls[r]])
        return out
    a, b = full(a0, a1), full(b0, b1)
    pa, pb = paf_cyc(a), paf_cyc(b)
    return pa[0] + pb[0] == 4 * w and all(pa[s] + pb[s] == 0 for s in range(1, 2 * w))

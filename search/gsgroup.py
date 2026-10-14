"""Finite abelian groups G = Z_{n1}^{k1} x ... with automorphism subgroups H, and the orbit tables used by gdfs.c.

Group spec: components separated by 'x', each "n" or "n^k" (Z_n^k); e.g. "25", "5^2", "3^2x11", "9^2".
Automorphism spec: one generator per '|', one part per component separated by '/', each part a unit "u" (k = 1) or a
k x k matrix "a,b;c,d" (rows separated by ';') acting on column vectors mod n.  Example for 3^2x11: "1,1;0,1/1|1,0;0,1/3".
Only automorphisms that act componentwise are supported; a matrix must be invertible mod n.

Elements are indexed in mixed radix over the coordinates.  For a block type (sym: one variable per pair of opposite
H-orbits with equal values, plus {0}; skew: pairs with opposite values, 0 excluded; plain: one variable per H-orbit)
and an x part U (a sequence over G, zero or H-invariant), `write_tables` writes the input of gdfs.c: per variable its
size, its PAF contributions at the representatives of the <H, -1>-orbits on G \\ {0}, the cross terms between
variables, and its contribution to X(s) = sum_g U(g) [W(g + s) + W(g - s)].
"""
import itertools
import struct

import numpy as np


def parse_group(spec):
    comps = []
    for part in spec.split('x'):
        if '^' in part:
            n, k = part.split('^')
            comps.append((int(n), int(k)))
        else:
            comps.append((int(part), 1))
    return comps


class Group:
    def __init__(self, spec):
        self.spec = spec
        self.comps = parse_group(spec)
        self.mods = [n for n, k in self.comps for _ in range(k)]
        self.t = int(np.prod(self.mods))
        self.d = len(self.mods)
        self.elts = np.array(list(itertools.product(*[range(m) for m in self.mods])), dtype=np.int64).reshape(self.t, self.d)
        self.mods_a = np.array(self.mods, dtype=np.int64)
        self.radix = np.array([int(np.prod(self.mods[i + 1:])) for i in range(self.d)], dtype=np.int64)
        # addition table and negation
        E = self.elts
        self.add = ((E[:, None, :] + E[None, :, :]) % self.mods_a) @ self.radix
        self.neg = ((-E) % self.mods_a) @ self.radix

    def index(self, coords):
        return int((np.asarray(coords) % self.mods_a) @ self.radix)

    def aut_perm(self, spec):
        """permutation (index -> index) of one automorphism given by its spec"""
        parts = spec.split('/')
        if len(parts) != len(self.comps):
            raise ValueError('automorphism spec does not match the components')
        img = np.zeros_like(self.elts)
        col = 0
        for (n, k), p in zip(self.comps, parts):
            if k == 1 and ';' not in p and ',' not in p:
                Mx = np.array([[int(p)]], dtype=np.int64)
            else:
                Mx = np.array([[int(x) for x in row.split(',')] for row in p.split(';')], dtype=np.int64)
            if Mx.shape != (k, k):
                raise ValueError('bad matrix shape')
            sub = self.elts[:, col:col + k]
            img[:, col:col + k] = (sub @ Mx.T) % n
            col += k
        perm = img @ self.radix
        if len(set(perm.tolist())) != self.t:
            raise ValueError(f'{spec} is not an automorphism')
        return perm

    def closure(self, perms):
        ident = tuple(range(self.t))
        H = {ident}
        fr = [ident]
        gens = [tuple(int(x) for x in p) for p in perms]
        while fr:
            x = fr.pop()
            for g in gens:
                y = tuple(g[i] for i in x)          # g o x
                if y not in H:
                    H.add(y)
                    fr.append(y)
        return [np.array(h, dtype=np.int64) for h in sorted(H)]


class Orbits:
    """H-orbits on G, <H,-1>-orbit representatives on G \\ {0} and the orbit-pair tables."""

    def __init__(self, G, H):
        self.G, self.H = G, H
        t = G.t
        if len(H) % 2 == 0:
            raise ValueError('H must have odd order')
        self.orb = -np.ones(t, dtype=np.int64)
        self.obs = []
        for g in range(t):
            if self.orb[g] < 0:
                o = sorted({int(h[g]) for h in H})
                for x in o:
                    self.orb[x] = len(self.obs)
                self.obs.append(o)
        self.no = len(self.obs)
        self.size = np.array([len(o) for o in self.obs], dtype=np.int64)
        self.negorb = np.array([self.orb[G.neg[o[0]]] for o in self.obs], dtype=np.int64)
        for i in range(1, self.no):
            if self.negorb[i] == i:
                raise ValueError('an H-orbit is closed under negation although |H| is odd')
        self.reps = []
        seen = set()
        for i in range(self.no):
            if self.obs[i][0] == 0 or i in seen:
                continue
            seen |= {i, int(self.negorb[i])}
            self.reps.append(self.obs[i][0])
        self.npr = len(self.reps)
        # Tp[o1, o2, k] = #{x in o1 : x + s_k in o2}, Tm with x - s_k
        Tp = np.zeros((self.no, self.no, self.npr), dtype=np.int64)
        Tm = np.zeros((self.no, self.no, self.npr), dtype=np.int64)
        for k, s in enumerate(self.reps):
            np.add.at(Tp, (self.orb, self.orb[G.add[:, s]], k), 1)
            np.add.at(Tm, (self.orb, self.orb[G.add[:, G.neg[s]]], k), 1)
        self.Tp, self.Tm = Tp, Tm

    def variables(self, btype, U=None, within=None):
        """variables of a block: list of [(orbit, coef), ...]; orbits carrying U are excluded; within: allowed orbits"""
        out, seen = [], set()
        for i in range(self.no):
            if i in seen:
                continue
            r0 = self.obs[i][0]
            if within is not None and not within[i]:
                continue
            if U is not None and U[r0]:
                continue
            if r0 == 0:
                if btype == 'skew':
                    continue
                out.append([(i, 1)])
                seen.add(i)
                continue
            j = int(self.negorb[i])
            if btype == 'plain':
                out.append([(i, 1)])
                seen.add(i)
            else:
                out.append([(i, 1), (j, 1 if btype == 'sym' else -1)])
                seen |= {i, j}
        return out

    def seq_from(self, var, vals):
        s = np.zeros(self.G.t, dtype=np.int64)
        for v, val in zip(var, vals):
            for o, c in v:
                s[self.obs[o]] = c * val
        return s

    def tables(self, var, mode, U=None):
        """(vsize, Dg, Cx, Lv) of the variables: vsize is the signed row-sum contribution (mode 0) or the weight
        (mode 1); Dg/Cx the PAF contributions at the representatives; Lv the contribution to X (needs U)."""
        nv, npr = len(var), self.npr
        if mode == 0:
            vsize = [int(sum(c * self.size[o] for o, c in v)) for v in var]
        else:
            vsize = [int(sum(self.size[o] for o, c in v)) for v in var]
        Dg = np.zeros((nv, npr), dtype=np.int64)
        Cx = np.zeros((nv, nv, npr), dtype=np.int64)
        Lv = np.zeros((nv, npr), dtype=np.int64)
        for j, v in enumerate(var):
            for o1, c1 in v:
                for o2, c2 in v:
                    Dg[j] += c1 * c2 * self.Tp[o1, o2]
            for i in range(j):
                for o1, c1 in v:
                    for o2, c2 in var[i]:
                        Cx[j, i] += c1 * c2 * (self.Tp[o1, o2] + self.Tp[o2, o1])
        if U is not None and np.any(U):
            # X(s_k) = sum_g U(g) [W(g + s_k) + W(g - s_k)]; per H-orbit o2 of W: L[o2, k]
            Uo = np.zeros(self.no, dtype=np.int64)
            for i in range(self.no):
                Uo[i] = U[self.obs[i][0]]
            L = np.einsum('a,abk->bk', Uo, self.Tp + self.Tm)
            for j, v in enumerate(var):
                for o, c in v:
                    Lv[j] += c * L[o]
        return vsize, Dg, Cx, Lv

    def write_tables(self, path, var, mode, target, scale=1, U=None):
        nv, npr = len(var), self.npr
        if nv > 64:
            raise ValueError(f'{nv} variables (max 64)')
        vsize, Dg, Cx, Lv = self.tables(var, mode, U)
        with open(path, 'wb') as f:
            f.write(struct.pack('5i', mode, nv, npr, target, scale))
            f.write(np.array(vsize, dtype=np.int32).tobytes())
            f.write(Dg.astype(np.int32).tobytes())
            f.write(Cx.astype(np.int32).tobytes())
            f.write(Lv.astype(np.int32).tobytes())
        return vsize

    def rep_weights(self):
        """number of shifts represented by each representative: the size of its <H,-1>-orbit"""
        return [2 * int(self.size[self.orb[s]]) for s in self.reps]

    def write_osa(self, path, blocks):
        """blocks: list of (mult, var, target, U) -> problem file for osa.c"""
        with open(path, 'wb') as f:
            f.write(struct.pack('2i', len(blocks), self.npr))
            for mult, var, target, U in blocks:
                vsize, Dg, Cx, Lv = self.tables(var, 0, U)
                f.write(struct.pack('3i', mult, len(var), target))
                f.write(np.array(vsize, dtype=np.int32).tobytes())
                f.write(Dg.astype(np.int32).tobytes())
                f.write(Cx.astype(np.int32).tobytes())
                f.write(Lv.astype(np.int32).tobytes())
            f.write(np.array(self.rep_weights(), dtype=np.int32).tobytes())


def read_list(path):
    with open(path, 'rb') as f:
        n, D = struct.unpack('2i', f.read(8))
        vec = np.frombuffer(f.read(2 * n * D), dtype=np.int16).reshape(n, D)
        rest = f.read()
    masks = np.frombuffer(rest, dtype=np.uint64)
    return n, D, vec, masks


def group_matrix(G, x):
    """X[g, h] = x(h - g)"""
    diff = G.add[:, G.neg].T          # diff[g, h] = h - g : add[h, neg[g]]
    return np.asarray(x)[diff]


def reflection(G):
    """R[g, h] = [g + h = 0]"""
    R = np.zeros((G.t, G.t), dtype=np.int64)
    R[np.arange(G.t), G.neg] = 1
    return R

"""Symmetry reduction of defect-position patterns (EXACT symmetries of the solution set of the K-level system): translation of B, C or D separately by a point of K
fixed by H_K, reflection of B, C or D separately, permutations of (B, C, D), units of K that normalise H_K acting on all four blocks.
Patterns are those of bbpat.enum_patterns: A = ((pair index, type), ...), B/C/D = ((orbit, type), ...).  Sym(sh).reduce(alive) keeps one pattern per class.
(Library module; used by run_pat.py, run_bb.py, val_sym.py.)"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
import itertools, numpy as np
from math import gcd


def unit_perms(sh):
    """list of permutations of the K-orbits (as maps orbit -> orbit) induced by the units of K (component-wise unit tuples; all normalise the abelian H_K), plus the map on A pairs"""
    Kg = sh.Kg; mods = Kg.mods
    ranges = [[a for a in range(1, m) if gcd(a, m) == 1] for m in mods]
    out = []; seen = set()
    sorb = sh.sorb; pos = sh.spos
    for u in itertools.product(*ranges):
        perm = Kg.mulperm(u)
        omap = tuple(int(sorb[int(perm[int(p[0])])]) for p in pos)
        # the unit must NORMALISE H_K (map orbits onto orbits): check every position of every orbit
        if any(len({int(sorb[int(perm[int(x)])]) for x in p}) != 1 for p in pos): continue
        if omap in seen: continue
        seen.add(omap); out.append(omap)
    return out


def fixed_points(sh):
    """points of K fixed by every generator of H_K"""
    N = sh.N; fp = np.ones(N, dtype=bool)
    for (u, ph) in sh.hgens:
        fp &= (np.asarray(u) == np.arange(N))
    return [int(x) for x in np.flatnonzero(fp)]


class Sym:
    def __init__(self, sh):
        self.sh = sh
        self.umaps = unit_perms(sh)
        self.fp = fixed_points(sh)
        add = sh.Kg.add
        pos = sh.spos; sorb = sh.sorb
        self.tmaps = []          # translation maps orbit -> orbit for g in fixed points
        for g in self.fp:
            self.tmaps.append(tuple(int(sorb[int(add[int(p[0]), g])]) for p in pos))
        self.negorb = tuple(int(sorb[int(sh.negperm[int(p[0])])]) for p in pos)
        # A pairs: pair j = (o, po); orbit map -> pair index and orientation
        self.pairof = {}
        for j, (o, po) in enumerate(sh.apairs): self.pairof[o] = (j, 0); self.pairof[po] = (j, 1)

    def canon_free(self, items):
        """canonical form of a free block pattern under translations by fixed points and reflection"""
        best = None
        for tm in self.tmaps:
            for refl in (False, True):
                img = []
                for (o, ty) in items:
                    o2 = tm[o]
                    if refl: o2 = self.negorb[o2]
                    img.append((o2, ty))
                img = tuple(sorted(img))
                if best is None or img < best: best = img
        return best

    def apply_unit(self, umap, A, B, C, D):
        # A: ((pair index, type)); a unit maps pair (o, po) to the pair containing umap[o]
        A2 = []
        for (j, ty) in A:
            o, po = self.sh.apairs[j]
            j2, side = self.pairof[umap[o]]
            A2.append((j2, ty))
        A2 = tuple(sorted(A2))
        free = [tuple((umap[o], ty) for (o, ty) in blk) for blk in (B, C, D)]
        return A2, free

    def canon(self, A, B, C, D):
        best = None
        for umap in self.umaps:
            A2, free = self.apply_unit(umap, A, B, C, D)
            cf = sorted(self.canon_free(b) for b in free)
            key = (A2, tuple(cf))
            if best is None or key < best: best = key
        return best

    def reduce(self, alive):
        seen = {}
        for a in alive:
            key = self.canon(a['A'], a['B'], a['C'], a['D'])
            if key not in seen: seen[key] = a
        return seen

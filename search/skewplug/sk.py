"""Skew Paley plug-in over G = K x F_q (q = 3 mod 4, a prime or the prime power 27), Goethals-Seidel: block A of skew type, B, C, D free.

Block X = s (x) I + t (x) Q + e (x) J on K (x0 = s + e on K x {0}, x1 = e + t on K x S, x2 = e - t on K x N; S = squares of F_q*).
K-level system (as for the real Goethals-Seidel plug-ins, gsk.py mode 'gs'):  (I) sum (s s~ + q t t~) = 4|G| delta,  (J) sum (q e e~ + s e~ + e s~ - t t~) = 0,  (Q) sum (t s~ - s t~) = 0   (q = 3 mod 4).
SKEW A:  A(-g) = -A(g) (g != 0), A(0) = 1.   -1 = non-square, so (k,f) -> (-k,-f) swaps the classes x1 <-> x2:
     x0(-k) = -x0(k) (k != 0), x0(0) = 1;   x2(-k) = -x1(k) (all k)
  =>  s(-k) = -s(k) (k != 0), s(0) = 1;  e(-k) = -e(k) (so e(0) = 0: NO defect at 0);  t(-k) = +t(k).
Defects of A come in pairs {k, -k} of the same type (p: s = 0; m: s = -2e).
"""
import numpy as np, itertools, math, sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from gsk import GShape, exact_paf_G, verify_g, enum_weight
from pk import KGrp, orbits_signed


class SShape(GShape):
    def __init__(self, Kmods, q, hunits=()):
        hg = [(tuple(u) if not np.isscalar(u) else u, 1) for u in hunits]
        super().__init__(Kmods, q, hg, 'gs')
        N = self.N
        assert N % 2 == 1 and q % 4 == 3
        self.hunits = [tuple(u) if not np.isscalar(u) else (int(u),) * len(self.Kmods) for u in hunits]
        orb, sig, forced, pos = orbits_signed(N, [(p, 1) for (p, ph) in self.hgens])
        self.sorb, self.spos, self.snu = orb, pos, len(pos)
        neg = self.negperm
        partner = [int(orb[neg[p[0]]]) for p in pos]
        zero = int(orb[0]); assert len(pos[zero]) == 1
        self.zero_orb = zero
        pairs = []; seen = set()
        for o in range(self.snu):
            if o == zero or o in seen: continue
            po = partner[o]; assert po != o, 'H_K contains -1'
            pairs.append((o, po)); seen |= {o, po}
        self.apairs = pairs
        self.apw = [len(pos[o]) for (o, po) in pairs]        # |O| (defect weight of a pair = 2|O|)
        # t-orbits of A: {0} and O u -O (sig +1)
        self.at_orbits = [[(0, 1)]] + [[(int(p), 1) for p in list(pos[o]) + list(pos[po])] for (o, po) in pairs]
        self.at_defect_idx = [None] + list(range(len(pairs)))
        self.nfree_units = self.kCD.nu

    # ----------------------------------------------------------------- A lists
    def enum_A(self, nA, mA=None, cap=60_000_000):
        """all skew-type (s,e) configs of A with defect weight nA (positions, = sum of 2|O| over defect pairs), optionally m-type weight mA.  Returns (S, E) int8 (M, N)."""
        N = self.N; pos = self.spos; pairs = self.apairs; np_ = len(pairs)
        w = [2 * len(pos[o]) for (o, po) in pairs]
        subsets = []
        def rec(j, rem, cur):
            if rem == 0: subsets.append(list(cur)); return
            if j == np_ or rem < 0: return
            if w[j] <= rem:
                cur.append(j); rec(j + 1, rem - w[j], cur); cur.pop()
            rec(j + 1, rem, cur)
        rec(0, nA, [])
        NOND = [(1, 0), (-1, 0)]
        DEF = [(0, 1), (-2, 1), (0, -1), (2, -1)]      # (s, e) on O ; on -O: (-s, -e)
        outS, outE = [], []; tot = 0
        for sub in subsets:
            ds = set(sub)
            opts = [DEF if j in ds else NOND for j in range(np_)]
            sizes = [len(o) for o in opts]
            M = int(np.prod(sizes, dtype=np.float64)); tot += M
            if tot > cap: raise RuntimeError('A list too large')
            idx = np.indices(sizes).reshape(np_, -1).T
            S = np.zeros((M, N), np.int8); E = np.zeros((M, N), np.int8)
            S[:, 0] = 1
            for j, (o, po) in enumerate(pairs):
                sel = np.array(opts[j], dtype=np.int8)[idx[:, j]]
                S[:, pos[o]] = sel[:, 0:1]; E[:, pos[o]] = sel[:, 1:2]
                S[:, pos[po]] = -sel[:, 0:1]; E[:, pos[po]] = -sel[:, 1:2]
            if mA is not None:
                mt = ((E != 0) & (S != 0)).sum(1)
                keep = mt == mA; S = S[keep]; E = E[keep]
            outS.append(S); outE.append(E)
        if not outS: return None
        S = np.concatenate(outS); E = np.concatenate(outE)
        return (S, E) if len(S) else None

    def A_weights(self):
        """achievable defect weights of A (subset sums of the pair weights 2|O|)"""
        reach = {0}
        for (o, po) in self.apairs:
            w = 2 * len(self.spos[o]); reach |= {x + w for x in reach}
        return sorted(reach)

    def A_dflags(self, e):
        """defect flags per A t-orbit for an e array (N,): zero orbit never defect; pair j defect iff e at the rep != 0"""
        out = [0]
        for (o, po) in self.apairs: out.append(int(e[int(self.spos[o][0])] != 0))
        return out

    def check_skew_A(self, s, t, e):
        neg = self.negperm; s = np.asarray(s); t = np.asarray(t); e = np.asarray(e)
        ok = s[0] == 1 and e[0] == 0
        k = np.arange(1, self.N)
        ok = ok and (s[neg[k]] == -s[k]).all() and (e[neg[k]] == -e[k]).all() and (t[neg] == t).all()
        return bool(ok)

    def check_solution(self, blocks):
        """exact K-level check: A skew in (s,t,e) form + (I),(J),(Q),(U) over all k"""
        sA, tA, eA = blocks[0]
        return self.check_skew_A(sA, tA, eA) and self.check_klevel(blocks)

    def g_target_skew_ok(self, seqs):
        """G-level exact check: A(0) = 1, A(-g) = -A(g), sum PAF = 4|G| delta"""
        mods = self.gmods(); n = self.G
        coords = np.array(np.unravel_index(np.arange(n), mods)).T
        m = np.array(mods); w = np.ones(len(mods), dtype=np.int64)
        for i in range(len(mods) - 2, -1, -1): w[i] = w[i + 1] * mods[i + 1]
        negi = (((-coords) % m) * w).sum(1)
        A = np.asarray(seqs[0], dtype=np.int64)
        okA = A[0] == 1 and all(A[negi[g]] == -A[g] for g in range(1, n))
        ok1 = all(set(np.unique(x).tolist()) <= {-1, 1} for x in seqs)
        paf = exact_paf_G(mods, seqs)
        tg = np.zeros(n, dtype=np.int64); tg[0] = 4 * n
        return bool(okA and ok1 and (paf == tg).all())


def block_full_from_se(sh, s, e, t):
    return (np.asarray(s, dtype=np.int64), np.asarray(t, dtype=np.int64), np.asarray(e, dtype=np.int64))

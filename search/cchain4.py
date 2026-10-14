#!/usr/bin/env python3
"""Exact compression-chain search for four-block Goethals-Seidel-type problems on Z_v with per-block symmetry:
blocks of kind "L" constant on the orbits of <H, -1> (symmetric), kind "A" on the orbits of H, and

    sum_k PAF_k(g) = T   for g != 0.

Modes:
  ssxx: v = N, T = 0, kinds L L A A (Szöllősi's complex Hadamard matrices of order 2N, search/ssxx.py);
  bsz:  v = N - 1, T = -4, kinds L L A A (bordered Szöllősi array, search/bsz.py).

Method (as search/chain.py for skew families). For d | v the compression x_d(j) = sum_{t = j mod d} x(t) of each block
is constant on the H-orbits of Z_d, and the compressions satisfy sum_k PAF_{x_k,d}(j) = T v/d (j != 0) on Z_d.
1. A dynamic program over the orbits lists, for each block, every joint compression (to Z_d1 and Z_d2) with its
   multiplicity; the Z_d1 compressions of the four blocks are joined exactly.
2. For each level-1 solution the Z_d2 compressions inside the four fibres are joined exactly.
3. For each consistent combination every member of the four fibres is listed (meet in the middle over the orbits)
   and the four lists are joined on all autocorrelation classes (search/join4x.c). Hits are rebuilt as complex
   Hadamard matrices and checked exactly.
Nothing is sampled; a run without truncation rules the structure/pattern out.

    python3 search/cchain4.py ssxx N g[,g2] [--pat a,b,c,d] [--d1 p --d2 q] [--budget S] [--log FILE]
    python3 search/cchain4.py bsz  N g[,g2] [--pat X] ...
"""
import argparse, itertools, json, math, os, struct, subprocess, sys, tempfile, time
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import workdir as W
import spence_lib as L
import chain as CH
import ssxx, bsz

KINDS = ["L", "L", "A", "A"]
MODES = {"ssxx": (["L", "L", "A", "A"], 0, False), "bsz": (["L", "L", "A", "A"], -4, True),
         "nw": (["A", "L", "L", "L"], 0, False), "bnw": (["A", "L", "L", "L"], -4, True)}


def closure(gens, v):
    S = {1}
    while True:
        T = S | {s * g % v for s in S for g in gens}
        if T == S:
            return sorted(S)
        S = T


class KChain:
    def __init__(self, v, gens, ds, T, kinds=KINDS):
        self.v, self.T, self.ds, self.kinds = v, T, list(ds), list(kinds)
        G = L.AbGroup([v]); self.G = G
        perms = [G.perm_from_matrix(np.array([[g]])) for g in gens]
        self.st = {"A": L.Structure(G, perms), "L": L.Structure(G, perms + [G.neg()])}
        assert np.array_equal(self.st["A"].cls, self.st["L"].cls)
        self.C = self.st["A"].C
        self.H = closure(gens, v)
        self.lev = []
        for d in self.ds:
            oid, reps, sizes = CH.h_orbits_mod(d, self.H)
            t = np.arange(v)
            Cm = {}
            for kind, st in self.st.items():
                M = np.zeros((st.O, len(reps)), dtype=np.int64)
                for r_i, r in enumerate(reps):
                    M[:, r_i] = np.bincount(st.orb[(t % d) == r], minlength=st.O)
                Cm[kind] = M
            self.lev.append((d, oid, reps, sizes, Cm))
        self.Ctot = {k: np.concatenate([l[4][k] for l in self.lev], axis=1) for k in self.st}
        self.splits = np.cumsum([0] + [len(l[2]) for l in self.lev])

    def states(self, kind, rowsum, maxstates):
        Ct = self.Ctot[kind].astype(np.int16)
        dim = Ct.shape[1]
        vt = np.dtype((np.void, 2 * dim))

        def dedup(K, cnt):
            key = np.ascontiguousarray(K).view(vt).ravel()
            u, idx, inv = np.unique(key, return_index=True, return_inverse=True)
            return K[idx], np.bincount(inv.ravel(), weights=cnt)

        K = (-Ct.sum(0))[None, :].astype(np.int16); cnt = np.ones(1)
        for o in range(self.st[kind].O):
            K, cnt = dedup(np.concatenate([K, K + 2 * Ct[o]]), np.concatenate([cnt, cnt]))
            if len(K) > maxstates:
                return None
        w = self.lev[0][3]
        rs = K[:, :len(w)].astype(np.int64) @ w
        keep = rs == rowsum
        return K[keep], cnt[keep]

    def paf(self, level, red):
        d, oid = self.lev[level][0], self.lev[level][1]
        X = red[:, oid]
        return np.stack([(X * np.roll(X, -k, axis=1)).sum(1) for k in range(1, d // 2 + 1)], axis=1)

    def target(self, level):
        d = self.lev[level][0]
        return np.full(d // 2, self.T * (self.v // d), dtype=np.int64)


class KHalves:
    """meet-in-the-middle tables of one block kind (all orbits free, base = all -1)"""

    def __init__(self, ch, kind):
        st, Ct = ch.st[kind], ch.Ctot[kind]
        self.items = list(range(st.O)); self.base = -Ct.sum(0); vecs = [2 * Ct[o] for o in self.items]
        m = len(self.items); h = m // 2; self.O = st.O
        self.R = np.random.default_rng(7).integers(1, 2 ** 62, size=Ct.shape[1], dtype=np.int64).astype(np.uint64)

        def half(lo, hi):
            n = hi - lo
            idx = np.arange(1 << n, dtype=np.int64)
            bits = ((idx[:, None] >> np.arange(n)) & 1).astype(np.int8)
            V = np.array(vecs[lo:hi], dtype=np.int64).reshape(n, -1)
            return bits, (bits.astype(np.int64) @ V).astype(np.int32)

        self.b1, self.s1 = half(0, h); self.b2, self.s2 = half(h, m)
        self.k1 = (self.s1.astype(np.int64).astype(np.uint64) * self.R).sum(1, dtype=np.uint64)
        k2 = (self.s2.astype(np.int64).astype(np.uint64) * self.R).sum(1, dtype=np.uint64)
        self.order = np.argsort(k2); self.ks = k2[self.order]

    def fibre(self, key, cap):
        need = np.asarray(key, dtype=np.int64) - self.base
        kn = (need.astype(np.uint64) * self.R).sum(dtype=np.uint64)
        want = kn - self.k1
        lo = np.searchsorted(self.ks, want, 'left'); hi = np.searchsorted(self.ks, want, 'right')
        cnt = hi - lo
        if int(cnt.sum()) > cap:
            return None
        parts = []
        nz = np.flatnonzero(cnt)
        for c0 in range(0, len(nz), 1 << 16):
            ii = nz[c0:c0 + (1 << 16)]
            k = cnt[ii]
            I = np.repeat(ii, k)
            off = np.arange(int(k.sum())) - np.repeat(np.cumsum(k) - k, k)
            J = self.order[lo[I] + off]
            ok = np.all(self.s1[I].astype(np.int64) + self.s2[J] == need, axis=1)
            I, J = I[ok], J[ok]
            if len(I) == 0:
                continue
            B = np.concatenate([self.b1[I], self.b2[J]], axis=1).astype(np.int8)
            parts.append((2 * B - 1).astype(np.int8))
        return np.concatenate(parts) if parts else np.zeros((0, self.O), dtype=np.int8)


def final_join(ch, lists, K=2):
    Vs = [np.rint(ch.st[kind].class_values(Y.astype(np.float64))).astype(np.int16) for kind, Y in zip(ch.kinds, lists)]
    sizes = [len(x) for x in Vs]
    perm = min(([0, 1, 2, 3], [0, 2, 1, 3], [0, 3, 1, 2]),
               key=lambda o: sizes[o[0]] * sizes[o[1]] + sizes[o[2]] * sizes[o[3]])
    fd, path = tempfile.mkstemp(suffix='.bin', dir=W.tmpdir()); os.close(fd)
    try:
        with open(path, 'wb') as f:
            f.write(struct.pack('5i', *[sizes[k] for k in perm], ch.C))
            f.write(np.full(ch.C, ch.T, dtype=np.int16).tobytes())
            for k in perm:
                f.write(np.ascontiguousarray(Vs[k]).tobytes())
        out = subprocess.run([W.bin('join4x'), path, str(K)], capture_output=True, text=True, check=True).stdout
    finally:
        os.remove(path)
    total = int(out.strip().split('TOTAL')[-1])
    sols = []
    for line in out.splitlines():
        if line.startswith('SOL'):
            ids = list(map(int, line.split()[1:]))
            q = [None] * 4
            for slot, k in enumerate(perm):
                q[k] = lists[k][ids[slot]]
            sols.append(q)
    return total, sols


def split_join(ch, lists, extra, K=2):
    """Final join through extra compression levels computed from the enumerated fibres: each list is grouped by its
    compressions to the extra levels, the compressed autocorrelations of the groups are joined exactly, and join4x
    runs only on consistent sub-lists. Exact."""
    if not extra:
        return final_join(ch, lists, K)
    groups = []
    for kind, Y in zip(ch.kinds, lists):
        keys = np.concatenate([Y.astype(np.int64) @ e.lev[0][4][kind] for e in extra], axis=1)
        u, inv = np.unique(keys, axis=0, return_inverse=True)
        groups.append((u, inv.ravel()))
    P = []
    for u, _ in groups:
        parts, a = [], 0
        for e in extra:
            h = e.lev[0][4]["A"].shape[1]
            parts.append(e.paf(0, u[:, a:a + h])); a += h
        P.append(np.concatenate(parts, axis=1))
    tgt = np.concatenate([e.target(0) for e in extra])
    try:
        quads = CH.join4(P, tgt)
    except CH.JoinTooLarge:
        return final_join(ch, lists, K)
    total, sols = 0, []
    for q in quads:
        sub = [lists[k][groups[k][1] == q[k]] for k in range(4)]
        t, s_ = final_join(ch, sub, K)
        total += t; sols += s_
        if sols:
            break
    return total, sols


def verify_sol(ch, mode, Ys):
    import bnw
    seqs = [Y.astype(np.int64)[ch.st[kind].orb] for kind, Y in zip(ch.kinds, Ys)]
    if mode == "ssxx":
        ssxx.complex_hadamard(seqs)
    elif mode == "bsz":
        bsz.build(seqs)
    elif mode == "nw":
        bnw.williamson_array(bnw.near_williamson_type(seqs)) if False else bnw.near_williamson_type(seqs)
    else:
        bnw.williamson_array(bnw.williamson_type(seqs))
    return seqs


def run(mode, N, gens, pat, d1=None, d2=None, maxfibre=3_000_000, maxstates=20_000_000, budget=3600, log=None,
        exhaustive=False, verbose=True):
    t0 = time.time()
    kinds, T, bordered = MODES[mode]
    v = N - 1 if bordered else N
    ps = [p for p in range(2, v) if v % p == 0 and all(p % r for r in range(2, int(p ** 0.5) + 1))]
    # prime-power divisors of v as candidate levels
    cand = sorted({p ** e for p in ps for e in range(1, 8) if v % (p ** e) == 0 and p ** e < v})
    if d1:
        levels = [d1] + ([d2] if d2 else [])
    else:
        counts1 = {}
        for d in cand:
            c1 = KChain(v, gens, [d], T, kinds)
            Sd = [c1.states(k, r, maxstates) for k, r in zip(kinds, pat)]
            if any(x is None for x in Sd):
                continue
            U = [np.unique(K_, axis=0).astype(np.int64) for K_, _ in Sd]
            try:
                n1 = 0 if any(len(u) == 0 for u in U) else len(CH.join4([c1.paf(0, u) for u in U], c1.target(0)))
            except CH.JoinTooLarge:
                continue
            counts1[d] = n1
            if n1 == 0:
                print(f"  RESULT {mode} N={N} H={gens} rows={pat}: ruled out at d={d} (empty level-1 join) "
                      f"[{time.time() - t0:.0f}s]", flush=True)
                return dict(status='ruled out', found=[])
        if not counts1:
            print("  state space too large at every level", flush=True); return dict(status='too-large', found=[])
        best = sorted(counts1, key=lambda p: counts1[p])
        levels = [best[0]]
        for d in best[1:]:
            if math.gcd(d, best[0]) == 1:
                levels.append(d); break
    ch = KChain(v, gens, levels, T, kinds)
    S = [ch.states(k, r, maxstates) for k, r in zip(kinds, pat)]
    if any(x is None for x in S):
        ch = KChain(v, gens, levels[:1], T, kinds); levels = levels[:1]
        S = [ch.states(k, r, maxstates) for k, r in zip(kinds, pat)]
        if any(x is None for x in S):
            print("  joint state space too large", flush=True); return dict(status='too-large', found=[])

    def project(level):
        a, b = ch.splits[level], ch.splits[level + 1]
        out = []
        for K_, c in S:
            u, inv = np.unique(K_[:, a:b], axis=0, return_inverse=True)
            out.append((u.astype(np.int64), inv.ravel(), np.bincount(inv.ravel(), weights=c)))
        return out

    proj = project(0)
    try:
        sols1 = CH.join4([ch.paf(0, u) for u, _, _ in proj], ch.target(0))
    except CH.JoinTooLarge:
        print("  level-1 join too large", flush=True); return dict(status='too-large', found=[])
    d1, d2 = levels[0], (levels[1] if len(levels) > 1 else 0)
    if verbose:
        print(f"{mode} N={N} H={gens} |H|={len(ch.H)} O_L={ch.st['L'].O} O_A={ch.st['A'].O} C={ch.C} rows={pat} "
              f"d1={d1} d2={d2}: states {[len(s[0]) for s in S]}, level-1 vectors {[len(p[0]) for p in proj]}, "
              f"level-1 solutions {len(sols1)} [{time.time() - t0:.0f}s]", flush=True)
    weight = [np.prod([proj[k][2][q[k]] for k in range(4)]) for q in sols1]
    branches = [sols1[i] for i in np.argsort(weight)[::-1]]
    halves = {}
    extra = [KChain(v, gens, [d], T, kinds) for d in cand if d not in levels and all(math.gcd(d, l) == 1 for l in levels)]
    # keep the extra levels coprime to each other as well (largest first)
    ex2 = []
    for e in sorted(extra, key=lambda e: -e.ds[0]):
        if all(math.gcd(e.ds[0], f.ds[0]) == 1 for f in ex2):
            ex2.append(e)
    extra = ex2
    truncated, found, combos, nsol, timed_out = 0, [], 0, 0, False
    if max(ch.st['A'].O, ch.st['L'].O) > 46:
        print('  too many orbits for the meet-in-the-middle expansion', flush=True); return dict(status='too-large', found=[])
    for bi, q in enumerate(branches):
        if time.time() - t0 > budget:
            truncated += 1; timed_out = True; break
        if d2:
            fib = []
            for k in range(4):
                rows = np.flatnonzero(proj[k][1] == q[k])
                fib.append((S[k][0][rows], S[k][1][rows]))
            a2, b2 = ch.splits[1], ch.splits[2]
            try:
                sols2 = CH.join4([ch.paf(1, Kf[:, a2:b2].astype(np.int64)) for Kf, _ in fib], ch.target(1))
            except CH.JoinTooLarge:
                truncated += 1; continue
            cands = [[fib[k][0][s[k]] for k in range(4)] for s in sols2]
            counts = [[fib[k][1][s[k]] for k in range(4)] for s in sols2]
            order = sorted(range(len(cands)), key=lambda i: -float(np.prod([float(c) for c in counts[i]])))
            cands = [cands[i] for i in order]; counts = [counts[i] for i in order]
        else:
            cands = [[S[k][0][np.flatnonzero(proj[k][1] == q[k])[0]] for k in range(4)]]
            counts = [[proj[k][2][q[k]] for k in range(4)]]
        for keys, cnt in zip(cands, counts):
            if time.time() - t0 > budget:
                truncated += 1; timed_out = True; break
            if max(cnt) > maxfibre:
                truncated += 1; continue
            lists = []
            for k, kind in enumerate(kinds):
                if kind not in halves:
                    halves[kind] = KHalves(ch, kind)
                Y = halves[kind].fibre(keys[k], maxfibre)
                if Y is None:
                    lists = None; break
                lists.append(Y)
            if lists is None:
                truncated += 1; continue
            combos += 1
            if any(len(Y) == 0 for Y in lists):
                continue
            total, sols = split_join(ch, lists, extra)
            nsol += total
            for Ys in sols:
                seqs = verify_sol(ch, mode, Ys)
                rec = {"N": N, "mode": mode, "gens": list(gens), "subgroup_order": len(ch.H),
                       "rowsums": [int(x.sum()) for x in seqs], "method": f"cchain4 d1={d1} d2={d2}",
                       "sequences": ["".join("+" if x > 0 else "-" for x in s) for s in seqs]}
                if mode == "bsz":
                    b = bsz.border(tuple(rec["rowsums"]))[0]
                    rec["border"] = [[int(z.real), int(z.imag)] for z in b]
                found.append(rec)
                if log:
                    with open(log, "a") as f:
                        f.write(json.dumps(rec) + "\n")
            if found and not exhaustive:
                break
        if verbose and (bi < 3 or bi % 50 == 0 or found):
            print(f"  branch {bi + 1}/{len(branches)}: combos {combos}, truncated {truncated}, solutions {nsol}, "
                  f"verified {len(found)} [{time.time() - t0:.0f}s]", flush=True)
        if found and not exhaustive:
            break
    status = 'found' if found else ('ruled out' if truncated == 0 else
                                    ('budget reached' if timed_out else f'nothing ({truncated} truncated)'))
    print(f"  RESULT {mode} N={N} H={gens} rows={pat}: {status}; branches {len(branches)}, combos {combos}, "
          f"solutions {nsol} [{time.time() - t0:.0f}s]", flush=True)
    return dict(status=status, found=found)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=list(MODES))
    ap.add_argument("N", type=int); ap.add_argument("gens")
    ap.add_argument("--pat"); ap.add_argument("--d1", type=int); ap.add_argument("--d2", type=int)
    ap.add_argument("--budget", type=float, default=3600); ap.add_argument("--log")
    ap.add_argument("--maxfibre", type=float, default=3e6); ap.add_argument("--maxstates", type=float, default=2e7)
    ap.add_argument("--exhaustive", action="store_true")
    a = ap.parse_args()
    gens = [int(x) for x in a.gens.split(",")]
    import bnw
    from bnw import Struct
    if a.mode == "ssxx":
        st = Struct(a.N, gens, plain=True)
        pats = [tuple(int(x) for x in a.pat.split(","))] if a.pat else ssxx.patterns(st)
    elif a.mode == "bsz":
        st = Struct(a.N, gens, plain=False)
        allp = bsz.patterns(st)
        pats = [p for n, p in allp if (not a.pat or n == a.pat)]
    else:
        st = Struct(a.N, gens, plain=(a.mode == "nw"))
        allp = bnw.patterns(st)
        pats = [p for n, p in allp if (not a.pat or n == a.pat)]
    for pat in pats:
        r = run(a.mode, a.N, gens, tuple(int(x) for x in pat), a.d1, a.d2, int(a.maxfibre), int(a.maxstates),
                a.budget, a.log, a.exhaustive)
        if r["found"] and not a.exhaustive:
            break


if __name__ == "__main__":
    main()

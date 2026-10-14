"""Symmetric Hadamard matrices of order 4v from propus-type quadruples (A, B, B, D) over non-cyclic abelian groups.

The array of propus.py (A symmetric, B and D arbitrary; or A and D symmetric) works over any abelian group G with
group-developed blocks X[g, h] = x(h - g) and R[g, h] = [g + h = 0]: every XR is symmetric, R X R = X^T, and the
blocks commute. The condition is PAF_A + 2 PAF_B + PAF_D = 4v delta on G with A(-g) = A(g).

For v with p^2 | v (gcd(p, v/p^2) = 1), G = Z_p x Z_p x ... x Z_m has automorphisms GL(2, p) on each split prime and
the units of Z_m on the rest (skew_group.py). The sweep takes, for every combination of a named subgroup of GL(2, p)
(Singer, diagonal, torus, unipotent, Borel, of any order) and a subgroup of Z_m^* (any order), the
structure H; A lives on the orbits of <H, -1>, B and D on the orbits of H (or of <H, -1>). Every row-sum pattern with
pools of at most MAXPOOL members is joined exactly with join4x, and every hit is rebuilt as a matrix and checked.

  python3 search/propus_group.py v [--mino 1] [--maxo 30] [--maxpool 3e7] [--maxwork W] [--K 1]
Hits go to $HADAMARD_WORK/sym_hits.jsonl with the group and the automorphism generators.
"""
import argparse, itertools, json, math, os, struct, subprocess, sys, tempfile, time
import numpy as np
from pathlib import Path as _P
sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
import spence_lib as L
import skew_group as SG
import skew_noncyclic as NC
import propus as PR
from gs_sweep_subgroups import subgroups


def unit_all(m):
    """All subgroups of Z_m^* (generators), including the trivial one."""
    if m == 1: return [("1", [])]
    out = [("1", [])]
    for H, gens in subgroups(m).items():
        if len(H) > 1: out.append((f"u{'.'.join(map(str, gens))}", list(gens)))
    return out


def candidates_all(p):
    """Subgroups of GL(2, p) of any order (a symmetric structure needs no odd order): Singer subgroups of every
    order d | p^2 - 1, diagonal subgroups <diag(x, x^k)>, the split torus, the unipotent group and Borel groups."""
    out = []
    S = NC.singer(p); N = p * p - 1
    for d in range(2, N + 1):
        if N % d == 0: out.append((f"singer{d}", [NC.mpow(S, N // d, p)]))
    g = next(x for x in range(2, p) if all(pow(x, (p - 1) // q, p) != 1 for q in range(2, p) if (p - 1) % q == 0 and all(q % r for r in range(2, q))))
    for d1 in [d for d in range(2, p) if (p - 1) % d == 0]:
        x = pow(g, (p - 1) // d1, p)
        for k in range(p - 1):
            out.append((f"diag{d1}_{k}", [np.array([[x, 0], [0, pow(x, k, p)]], dtype=np.int64)]))
    out.append(("torus", [np.array([[g, 0], [0, 1]]), np.array([[1, 0], [0, g]])]))
    U = np.array([[1, 1], [0, 1]], dtype=np.int64)
    out.append(("unip", [U]))
    for d in [d for d in range(2, p) if (p - 1) % d == 0]:
        x = pow(g, (p - 1) // d, p)
        for k in range(d):
            out.append((f"borel{d}_{k}", [U, np.array([[x, 0], [0, pow(x, k, p)]], dtype=np.int64)]))
    return out


def build(v, ps, m, choice):
    """Group and automorphism matrices (full k x k) of the direct product of the chosen parts."""
    mods = []
    for p in ps: mods += [p, p]
    if m > 1: mods.append(m)
    G = L.AbGroup(mods)
    k = len(mods); mats = []; off = 0
    for p, (name, gens) in zip(ps, choice[:len(ps)]):
        for Mx in gens:
            A = np.eye(k, dtype=np.int64); A[off:off + 2, off:off + 2] = np.asarray(Mx) % p; mats.append(A)
        off += 2
    if m > 1:
        for u in choice[-1][1]:
            A = np.eye(k, dtype=np.int64); A[off, off] = u % m; mats.append(A)
    return G, mats


def structures(v, mino, maxo):
    ps, m = SG.split(v)
    if not ps: return
    parts = [[("I", [])] + candidates_all(p) for p in ps] + [unit_all(m)]
    seen = set()
    for choice in itertools.product(*parts):
        if all(not c[1] for c in choice): continue
        G, mats = build(v, ps, m, list(choice))
        try:
            perms = [G.perm_from_matrix(A) for A in mats]
        except ValueError:
            continue
        st_b = L.Structure(G, perms)
        st_ad = L.Structure(G, perms + [G.neg()])
        if not (mino <= max(st_b.O, st_ad.O) <= maxo): continue
        key = (tuple(sorted(np.bincount(st_b.orb).tolist())), tuple(sorted(np.bincount(st_ad.orb).tolist())), st_b.C)
        if key in seen: continue
        seen.add(key)
        yield "+".join(c[0] for c in choice), G, mats, st_ad, st_b


def group_paf(G, x):
    """PAF on G: P(g) = sum_h x(h) x(h + g), for every g (index order)."""
    return np.array([int(np.dot(x, x[G.translate(g)])) for g in range(G.v)], dtype=np.int64)


def propus_matrix(G, A, B, D):
    v = G.v
    diff = G.index(G.coords[None, :, :] - G.coords[:, None, :])        # diff[g, h] = h - g
    cA, cB, cD = A[diff], B[diff], D[diff]
    R = np.zeros((v, v), dtype=np.int64); R[np.arange(v), G.neg()] = 1
    BR, DR = cB @ R, cD @ R
    return np.block([[cA, cB, cB, cD], [cB.T, DR, -cA, -BR], [cB.T, -cA, -DR, BR], [cD.T, -BR, BR, -cA]])


def search(v, name, G, mats, st_ad, st_b, symB, freeD, maxpool, K, log, maxwork=float('inf')):
    st_a = st_ad
    st_bb = st_ad if symB else st_b
    st_d = st_b if freeD else st_ad
    pm = v % 4
    pats = PR.patterns(v)
    log(f"v={v} {name} symB={symB} freeD={freeD}: orbits A {st_a.O}, B {st_bb.O}, D {st_d.O}, classes {st_ad.C}, patterns {len(pats)}")
    found = []
    for ra, rb, rd in pats:
        if not freeD and ra > rd: continue
        sizes = [st_a.pool(ra), st_bb.pool(rb), st_d.pool(rd)]
        if min(sizes) == 0 or max(sizes) > maxpool:
            log(f"  rows ({ra}, {rb}, {rd}): pools {[int(s) for s in sizes]} skipped"); continue
        work = float(np.prod(sorted(sizes)[:2]))                # pairs of the two smallest lists, probed against the third
        if work > maxwork:
            log(f"  rows ({ra}, {rb}, {rd}): pools {[int(s) for s in sizes]} work {work:.1e} skipped"); continue
        t = time.time()
        def vals(st, Y, scale, chunk=1 << 18):
            out = np.empty((len(Y), st.C), dtype=np.int16)
            for c0 in range(0, len(Y), chunk):
                out[c0:c0 + chunk] = scale * np.rint((st.class_values(Y[c0:c0 + chunk].astype(np.float64)) - pm) / 4)
            return out
        Ya = PR.pool_vectors(st_a, ra); Yb = PR.pool_vectors(st_bb, rb)
        Yd = Ya if (rd == ra and not freeD) else PR.pool_vectors(st_d, rd)
        Va = vals(st_a, Ya, 1); Vb = vals(st_bb, Yb, 2); Vd = Va if (rd == ra and not freeD) else vals(st_d, Yd, 1)
        Z = np.zeros((1, st_ad.C), dtype=np.int16)
        roles = sorted([(len(Va), 0), (len(Vd), 1), (len(Vb), 2)])   # the two smallest lists form the pairs
        byrole = {0: Va, 1: Vd, 2: Vb}
        lists = [byrole[roles[0][1]], byrole[roles[1][1]], byrole[roles[2][1]], Z]
        slot = {r: i for i, (_, r) in enumerate(roles)}            # role (0 = A, 1 = D, 2 = B) -> position in lists
        fd, path = tempfile.mkstemp(suffix='.bin', dir=W.tmpdir()); os.close(fd)
        try:
            with open(path, 'wb') as f:
                f.write(struct.pack('iiiii', *[len(x) for x in lists], st_ad.C))
                f.write(np.full(st_ad.C, -pm, dtype=np.int16).tobytes())
                for x in lists: f.write(np.ascontiguousarray(x).tobytes())
            out = subprocess.run([W.bin('join4x'), path, str(K)], capture_output=True, text=True, check=True).stdout
        finally:
            os.remove(path)
        total = int(out.strip().split('TOTAL')[-1])
        log(f"  rows ({ra}, {rb}, {rd}): pools {[int(s) for s in sizes]}: {total} solutions [{time.time() - t:.1f}s]")
        for line in out.splitlines():
            if not line.startswith('SOL'): continue
            ids = list(map(int, line.split()[1:]))
            ia, idd, ib = ids[slot[0]], ids[slot[1]], ids[slot[2]]
            A = Ya[ia].astype(np.int64)[st_a.orb]; B = Yb[ib].astype(np.int64)[st_bb.orb]; D = Yd[idd].astype(np.int64)[st_d.orb]
            tot = group_paf(G, A) + 2 * group_paf(G, B) + group_paf(G, D)
            if not (tot[0] == 4 * v and not tot[1:].any()): log('  WRONG autocorrelations'); continue
            if not np.array_equal(A, A[G.neg()]): log('  A not symmetric'); continue
            Hm = propus_matrix(G, A, B, D)
            ok = np.array_equal(Hm, Hm.T) and np.array_equal(Hm @ Hm.T, 4 * v * np.eye(4 * v, dtype=Hm.dtype))
            log(f"  hit: symmetric Hadamard {ok}")
            if ok:
                rec = dict(n=v, group=[int(x) for x in G.mods], structure=name, automorphisms=[A_.tolist() for A_ in mats],
                           symB=symB, freeD=freeD, rows=[int(A.sum()), int(B.sum()), int(D.sum())],
                           seqs=[''.join('+' if x > 0 else '-' for x in s) for s in (A, B, D)])
                with open(W.hits('sym_hits.jsonl'), 'a') as fh: fh.write(json.dumps(rec) + '\n')
                found.append(rec)
        if found: break
    return found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('v', type=int); ap.add_argument('--mino', type=int, default=1); ap.add_argument('--maxo', type=int, default=30)
    ap.add_argument('--maxpool', type=float, default=3e7); ap.add_argument('--K', type=int, default=1)
    ap.add_argument('--maxwork', type=float, default=float('inf'), help='skip row patterns whose pair work exceeds this')
    ap.add_argument('--only-freeD', action='store_true', help='sweep only the variant with A symmetric and B, D free')
    a = ap.parse_args()
    log = lambda s: print(s, flush=True)
    jobs = sorted(structures(a.v, a.mino, a.maxo), key=lambda j: j[3].O + j[4].O)
    log(f"v={a.v}: {len(jobs)} non-cyclic structures")
    for name, G, mats, st_ad, st_b in jobs:
        for symB, freeD in (((False, True),) if a.only_freeD else ((True, False), (False, False), (False, True))):
            if search(a.v, name, G, mats, st_ad, st_b, symB, freeD, a.maxpool, a.K, log, a.maxwork):
                log(f"HIT v={a.v} {name} symB={symB} freeD={freeD}"); sys.exit(0)


if __name__ == '__main__':
    main()

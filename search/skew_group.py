"""Skew Goethals-Seidel families over non-cyclic abelian groups G = (Z_p x Z_p)^... x Z_m.

For v with p^2 | v, Z_v can be replaced by G = Z_p x Z_p x Z_{v/p^2} (gcd(p, v/p^2) = 1). Automorphisms
are block-diagonal: an element of GL(2, p) on each split prime and a unit of Z_m on the rest. Every H of
odd order admits a skew first block, since no element of odd order has eigenvalue -1. Group-developed
blocks A(g, h) = a(h - g) with R(g, h) = [g + h = 0] give a Goethals-Seidel array over any abelian group.

  python3 skew_group.py rank v [MINO MAXO]             # structures, Gaussian E, sampled yield, exact work
  python3 skew_group.py run v NAME b,c,d [CAP]         # exact search (find4b) of one structure
"""
import itertools, json, math, os, struct, subprocess, sys, tempfile, time
import numpy as np
from pathlib import Path as _P; sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
import spence_lib as L
from gs_sweep_subgroups import subgroups
import skew_rank as R
from skew_ktree import pool_exact
import skew_noncyclic as NC
HERE = str(W.HERE)


def odd(n):
    while n % 2 == 0: n //= 2
    return n


def split(v):
    """Split primes p with p^2 | v and gcd(p, v/p^2) = 1; returns (list of p, m)."""
    ps, m = [], v
    for p in range(3, int(math.isqrt(v)) + 1, 2):
        if all(p % q for q in range(3, int(p ** 0.5) + 1, 2)) and m % (p * p) == 0 and (m // (p * p)) % p != 0:
            ps.append(p); m //= p * p
    return ps, m


def gl2_odd(p):
    """Named odd-order generating sets in GL(2, p) (lists of 2x2 matrices)."""
    out = [("I", [])]
    for name, gens in NC.candidates(p):
        out.append((name, gens))
    return out


def unit_odd(m):
    """Odd-order subgroups of Z_m^* (generators), including the trivial one."""
    if m == 1: return [("1", [])]
    out = [("1", [])]
    for H, gens in subgroups(m).items():
        if len(H) % 2 == 1 and len(H) > 1: out.append((f"u{'.'.join(map(str, gens))}", list(gens)))
    return out


def build(v, ps, m, choice):
    """choice: list of (name, gens) per split prime, then (name, gens) for Z_m. Direct product of the parts."""
    mods = []
    for p in ps: mods += [p, p]
    if m > 1: mods.append(m)
    G = L.AbGroup(mods)
    k = len(mods)
    perms = []
    off = 0
    for p, (name, gens) in zip(ps, choice[:len(ps)]):
        for Mx in gens:
            A = np.eye(k, dtype=np.int64); A[off:off + 2, off:off + 2] = np.asarray(Mx) % p
            perms.append(G.perm_from_matrix(A))
        off += 2
    if m > 1:
        for u in choice[-1][1]:
            A = np.eye(k, dtype=np.int64); A[off, off] = u % m
            perms.append(G.perm_from_matrix(A))
    if not perms: perms = [np.arange(G.v)]
    return G, L.Structure(G, perms)


def structures(v, mino, maxo):
    ps, m = split(v)
    if not ps: return
    parts = [gl2_odd(p) for p in ps] + ([unit_odd(m)] if m > 1 else [])
    seen = set()
    for choice in itertools.product(*parts):
        if all(not c[1] for c in choice): continue
        try:
            G, st = build(v, ps, m, list(choice) if m > 1 else list(choice) + [("1", [])])
        except ValueError:
            continue
        if not (mino <= st.O <= maxo): continue
        key = (tuple(sorted(np.bincount(st.orb).tolist())), st.C, tuple(sorted(np.bincount(st.cls[st.cls >= 0]).tolist())))
        if key in seen: continue
        seen.add(key)
        yield "+".join(c[0] for c in choice), G, st


def gs_group(seqs, G):
    """Goethals-Seidel array over the abelian group G from four +-1 functions on G (index order)."""
    v = G.v
    diff = G.index(G.coords[None, :, :] - G.coords[:, None, :])        # diff[g, h] = h - g
    Ms = [s[diff] for s in seqs]
    Rm = np.zeros((v, v), dtype=np.int64); Rm[np.arange(v), G.neg()] = 1
    A, B, C, D = Ms
    return np.block([[A, B @ Rm, C @ Rm, D @ Rm],
                     [-B @ Rm, A, D.T @ Rm, -C.T @ Rm],
                     [-C @ Rm, -D.T @ Rm, A, B.T @ Rm],
                     [-D @ Rm, C.T @ Rm, -B.T @ Rm, A]])


def patterns(st, v):
    lim = int(math.isqrt(4 * v)); out = []
    for b in range(1, lim + 1, 2):
        for c in range(b, lim + 1, 2):
            d2 = 4 * v - 1 - b * b - c * c
            if d2 < c * c: continue
            d = math.isqrt(d2)
            if d * d == d2 and all(st.pool(r) > 0 for r in (b, c, d)): out.append((b, c, d))
    return out


def exact(v, G, st, q, cap=6e10, K=2):
    pm = v % 4
    rep = np.array([np.flatnonzero(st.orb == o)[0] for o in range(st.O)])
    partner = st.orb[G.neg()[rep]]
    pairs = [(o, partner[o]) for o in range(1, st.O) if o < partner[o]]
    s = 1 - 2 * ((np.arange(2 ** len(pairs))[:, None] >> np.arange(len(pairs))) & 1)
    YA = np.zeros((len(s), st.O)); YA[:, 0] = 1
    for j, (o, p_) in enumerate(pairs): YA[:, o] = s[:, j]; YA[:, p_] = -s[:, j]
    lists = [YA] + [pool_exact(st, r).astype(np.float64) for r in q]
    sizes = [len(x) for x in lists]
    order = [0] + sorted(range(1, 4), key=lambda k: sizes[k])        # A with the largest; the other two right
    order = [0, order[1], order[2], order[3]]
    work = sizes[0] * sizes[order[3]] + sizes[order[1]] * sizes[order[2]]
    if work > cap: return None, work
    Vs = []
    for k in order:
        Y = lists[k]; V = np.empty((len(Y), st.C), np.int16)
        step = max(1000, int(1.5e8 // (st.C * st.O)))
        for s0 in range(0, len(Y), step): V[s0:s0 + step] = np.rint((st.class_values(Y[s0:s0 + step]) - pm) / 4)
        Vs.append(V)
    fd, path = tempfile.mkstemp(suffix='.bin', dir=W.tmpdir()); os.close(fd)
    try:
        with open(path, 'wb') as f:
            f.write(struct.pack('iiiii', *[len(V) for V in Vs], st.C)); f.write(np.full(st.C, -pm, dtype=np.int16).tobytes())
            for V in Vs: f.write(np.ascontiguousarray(V).tobytes())
        out = subprocess.run([W.bin('find4b'), path, str(K)], capture_output=True, text=True, check=True).stdout
    finally:
        os.remove(path)
    total = int(out.strip().split('TOTAL')[-1])
    found = []
    for line in out.splitlines():
        if not line.startswith('SOL'): continue
        ids = list(map(int, line.split()[1:]))
        seqs = [lists[order[k]][ids[k]].astype(np.int64)[st.orb] for k in range(4)]
        Hm = gs_group(seqs, G).astype(np.float64)
        had = np.array_equal(Hm @ Hm.T, 4 * v * np.eye(4 * v)); skew = np.array_equal(Hm + Hm.T, 2 * np.eye(4 * v))
        found.append((had, skew, seqs))
    return (total, found), work


if __name__ == '__main__' and sys.argv[1] != 'certify':
    v = int(sys.argv[2])
    if sys.argv[1] == 'rank':
        mino, maxo = (int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else (10, 60)
        rng = np.random.default_rng(1)
        for name, G, st in structures(v, mino, maxo):
            res = R.rank_structure(st, v, 1.2e7, rng=rng)
            for lE2, lE, q, lN in sorted(res, reverse=True)[:2]:
                NA = 10 ** lN[0]; Ns = sorted(10 ** x for x in lN[1:])
                work = NA * Ns[2] + Ns[0] * Ns[1]
                print(f"{lE2:7.2f} v={v} mods={list(G.mods)} {name} O={st.O} C={st.C} rows=(1,{q}) log10E={lE:.2f} "
                      f"log10N=({','.join(f'{x:.1f}' for x in lN)}) work={work:.1e}", flush=True)
    elif sys.argv[1] == 'ktree':
        import skew_ktree as K
        name = sys.argv[3]; q = [int(x) for x in sys.argv[4].split(',')]
        maxtau = int(sys.argv[5]) if len(sys.argv) > 5 else 2000
        for nm, G, st in structures(v, 1, 10 ** 9):
            if nm != name: continue
            found = K.run(v, [], maxtau=maxtau, top=8, only=q, sols=2, maxpool=1.2e7, prebuilt=(G, st))
            for h in found:
                h.update(mods=[int(x) for x in G.mods], structure=name)
                with open(W.hits('skew_group_hits.jsonl'), 'a') as f: f.write(json.dumps(h) + '\n')
            print(f"done: {len(found)} verified", flush=True)
            break
    else:
        name, q = sys.argv[3], tuple(int(x) for x in sys.argv[4].split(','))
        cap = float(sys.argv[5]) if len(sys.argv) > 5 else 6e10
        for nm, G, st in structures(v, 1, 10 ** 9):
            if nm != name: continue
            t = time.time()
            res, work = exact(v, G, st, q, cap)
            print(f"v={v} {name} mods={list(G.mods)} O={st.O} C={st.C} rows=(1,{q}) work={work:.1e}: "
                  f"{'skipped' if res is None else str(res[0]) + ' solutions'} [{time.time() - t:.0f}s]", flush=True)
            if res:
                for had, skew, seqs in res[1]:
                    print(f"  VERIFIED Hadamard={had} skew={skew}", flush=True)
                    if had and skew:
                        with open(W.hits('skew_group_hits.jsonl'), 'a') as f:
                            f.write(json.dumps(dict(v=v, mods=[int(x) for x in G.mods], structure=name, rows=[1, *q],
                                                    seqs=[''.join('+' if x > 0 else '-' for x in s) for s in seqs])) + '\n')
            break


def automorphism_matrices(v, name):
    """Integer matrices (on all coordinates) of the generators of the named structure."""
    ps, m = split(v)
    parts = [dict(gl2_odd(p)) for p in ps] + ([dict(unit_odd(m))] if m > 1 else [])
    names = name.split('+')
    k = 2 * len(ps) + (1 if m > 1 else 0)
    mats, off = [], 0
    for p, nm, part in zip(ps, names, parts):
        for Mx in part[nm]:
            A = np.eye(k, dtype=np.int64); A[off:off + 2, off:off + 2] = np.asarray(Mx) % p; mats.append(A)
        off += 2
    if m > 1:
        for u in parts[-1][names[-1]]:
            A = np.eye(k, dtype=np.int64); A[off, off] = u % m; mats.append(A)
    return mats


def certify(v, k=0):
    import gzip, hashlib
    sys.path.insert(0, str(W.ROOT / 'search'))
    hits = [json.loads(l) for l in open(W.hits('skew_group_hits.jsonl'))]
    h = [x for x in hits if x['n'] == v and x['hadamard'] and x['skew']][k]
    mods = h['mods']
    mats = [np.array(h['generator'], dtype=np.int64)] if 'generator' in h else automorphism_matrices(v, h['structure'])
    G = L.AbGroup(mods)
    # subgroup order: closure of the generator permutations
    perms = [G.perm_from_matrix(M) for M in mats]
    seen = {tuple(range(G.v))}; frontier = [np.arange(G.v)]
    while frontier:
        q = frontier.pop()
        for p in perms:
            r = tuple(p[q])
            if r not in seen: seen.add(r); frontier.append(np.array(r))
    seqs = [np.array([1 if c == '+' else -1 for c in s], dtype=np.int64) for s in h['seqs']]
    Hm = gs_group(seqs, G)
    import io
    buf = io.BytesIO(); np.savetxt(buf, Hm, fmt="%d", delimiter=","); data = buf.getvalue()
    order = 4 * v; out = os.path.join(str(W.ROOT / 'skew'), str(order))
    if os.path.exists(os.path.join(out, 'certificate.json')):
        sys.exit(f"{out}/certificate.json exists; not overwritten")
    os.makedirs(out, exist_ok=True)
    cert = {"order": order, "construction": "goethals-seidel-group", "skew": True, "group": [int(x) for x in mods],
            "structure": h['structure'], "automorphisms": [M.tolist() for M in mats], "subgroup_order": len(seen),
            "rowsums": [int(x.sum()) for x in seqs], "sequences": h['seqs'], "sha256": hashlib.sha256(data).hexdigest()}
    with open(os.path.join(out, 'certificate.json'), 'w') as f: json.dump(cert, f, indent=1)
    with open(os.path.join(out, 'sequences.txt'), 'w') as f: f.write(''.join(s + '\n' for s in h['seqs']))
    with gzip.open(os.path.join(out, f'H{order}.csv.gz'), 'wb', compresslevel=9) as f: f.write(data)
    import verify
    from pathlib import Path
    n_, digest = verify.check(Path(out))
    print(f"skew order {n_}: group {mods}, |H| = {len(seen)}, rows {cert['rowsums']}, sha256 {digest}")


if __name__ == '__main__' and sys.argv[1] == 'certify':
    certify(int(sys.argv[2]))

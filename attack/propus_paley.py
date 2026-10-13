"""Symmetric Hadamard matrices from propus quadruples on symmetric Paley products, by configurations over Z_m.

Structure: v = m q with q = 1 (mod 4) a prime (or a prime power p^k: the group is then Z_m x GF(q)), gcd(m, q) = 1,
H = {h = 1 (mod m), h a square in GF(q)}. Since -1 is a square, chi (the quadratic character of GF(q)) is symmetric.
Every H-invariant block is X = u (x) delta + w (x) J + z (x) chi with u, w, z on Z_m, and at each position either
w = 0 and u, z = +-1 (z-type) or z = 0, w = +-1 and u in {0, -2w} (w-type). With chi chi = q delta - J, J chi = 0,
J J = q J, the propus condition sum_k c_k X_k X_k^(-1) = 4mq delta (c = 1, 2, 1 for A, B, D; A symmetric) becomes
    (I)   sum_k c_k [u_k u_k* + q z_k z_k*]                = 4mq delta
    (II)  sum_k c_k [u_k w_k* + w_k u_k* + q w_k w_k*]    = sum_k c_k z_k z_k*
    (III) sum_k c_k [u_k z_k* + z_k u_k*]                 = 0
over Z[Z_m]. The identity coefficient of (II) gives (q + 1) N_w - 4 N_- = 4m for the c-weighted numbers of w-positions
(N_w) and of those with u = -2w (N_-), so N_w is even. With N_w = 2 the w-positions are
    B: one w-position (translated to 0): q = 2m - 1 (u = 0) or q = 2m + 3 (u = -2w);
    A and D: one each, at 0 (A symmetric): q = 2m - 1, 2m + 1 or 2m + 3 (N_- = 0, 1, 2).
(Two w-positions in one block leave a term of size q^2 in the autocorrelation of the compressions.)

The search, exact and without decoupling assumptions:
  stage 1: the compressions y_k = u_k + q w_k with PAF_yA + 2 PAF_yB + PAF_yD = 0 off 0 (y_A symmetric);
  stage 2: for each stage-1 triple, z_A (symmetric), z_B, z_D (zero at the w-positions) with
           PAF_zA + 2 PAF_zB + PAF_zD = R and the c-weighted cross terms (III) = 0, R(s) = sum_k c_k w_k(0)(u_k(s) + u_k(-s)).
Both stages are exact joins (join4x). Hits are lifted to Z_m x GF(q), checked as symmetric Hadamard matrices and
recorded in $HADAMARD_WORK/sym_hits.jsonl.

  python3 attack/propus_paley.py v [--q Q] [--field p,k] [--maxy N] [--nsol N]
"""
import argparse, json, math, sys, time
import numpy as np
from pathlib import Path as _P
sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
import verify
import paley_config as PC
import spence_lib as L
import propus_group as PG


def pafs0(X, h):
    """PAF(s), s = 1..h."""
    return PC.pafs(X, h)


def cross(Z, u, h):
    """sum_a [u(a) z(a+s) + z(a) u(a+s)], s = 0..h (symmetric in s)."""
    return np.stack([(np.roll(Z, -s, axis=1) * u[None, :] + Z * np.roll(u, -s)[None, :]).sum(1) for s in range(0, h + 1)], axis=1)


def symmetric_vecs(m, at0=None):
    """All symmetric +-1 sequences on Z_m (value at 0 fixed to at0 if given)."""
    h = (m - 1) // 2
    S = PC.allpm(h + (0 if at0 is not None else 1))
    X = np.zeros((len(S), m), np.int64)
    if at0 is None:
        X[:, 0] = S[:, 0]; T = S[:, 1:]
    else:
        X[:, 0] = at0; T = S
    for i in range(1, h + 1): X[:, i] = T[:, i - 1]; X[:, m - i] = T[:, i - 1]
    return X


def free_vecs(m, at0=None, zero_at=()):
    """All +-1 sequences on Z_m; position 0 fixed to at0 if given; positions in zero_at set to 0 (z-lists)."""
    if at0 is not None or zero_at:
        S = PC.allpm(m - 1); X = np.zeros((len(S), m), np.int64); X[:, 1:] = S
        X[:, 0] = 0 if at0 is None else at0
        return X
    return PC.allpm(m).astype(np.int64)


def layouts(m, q):
    """(name, wA, wB, wD, uA0, uB0, uD0) for N_w = 2: w = 1 at position 0 of the named blocks, u there 0 or -2."""
    out = []
    if q == 2 * m - 1: out += [("B", 0, 1, 0, None, 0, None), ("AD", 1, 0, 1, 0, None, 0)]
    if q == 2 * m + 3: out += [("B", 0, 1, 0, None, -2, None), ("AD", 1, 0, 1, -2, None, -2)]
    if q == 2 * m + 1: out += [("AD", 1, 0, 1, -2, None, 0), ("AD", 1, 0, 1, 0, None, -2)]
    return out


def chi_table(q, field):
    if field:
        chi, M = PC.gf_model(*field)
        return chi, M
    chi = np.zeros(q, np.int64)
    for x in range(1, q): chi[x * x % q] = 1
    chi[(chi == 0) & (np.arange(q) != 0)] = -1
    return chi, None


def lift(m, q, field, U, Wv, Z):
    chi, M = chi_table(q, field)
    if field:      # Z_m x GF(q) in mixed radix (Z_m slowest)
        a = np.repeat(np.arange(m), q); b = np.tile(np.arange(q), m)
        return [(u[a] * (b == 0) + w[a] + z[a] * chi[b]).astype(np.int64) for u, w, z in zip(U, Wv, Z)], M
    v = m * q; t = np.arange(v); a, b = t % m, t % q
    return [(u[a] * (b == 0) + w[a] + z[a] * chi[b]).astype(np.int64) for u, w, z in zip(U, Wv, Z)], None


def check(m, q, field, seqs):
    A, B, D = seqs
    if field:
        G = L.AbGroup([m] + [field[0]] * field[1])
        tot = PG.group_paf(G, A) + 2 * PG.group_paf(G, B) + PG.group_paf(G, D)
        if not (tot[0] == 4 * m * q and not tot[1:].any()) or not np.array_equal(A, A[G.neg()]): return False
        Hm = PG.propus_matrix(G, A, B, D)
    else:
        v = m * q
        tot = verify.paf_sum([A, B, D]) + verify.paf_sum([B])
        if not (tot[0] == 4 * v and not tot[1:].any()) or not np.array_equal(A, A[(-np.arange(v)) % v]): return False
        cA, cB, cD = (verify.circulant(x) for x in (A, B, D))
        R = verify.back_identity(v); BR, DR = cB @ R, cD @ R
        Hm = np.block([[cA, cB, cB, cD], [cB.T, DR, -cA, -BR], [cB.T, -cA, -DR, BR], [cD.T, -BR, BR, -cA]])
    n = Hm.shape[0]
    return np.array_equal(Hm, Hm.T) and np.array_equal(Hm @ Hm.T, n * np.eye(n, dtype=Hm.dtype))


def search(v, m, q, field, lay, a, log):
    name, wA, wB, wD, uA0, uB0, uD0 = lay
    h = (m - 1) // 2; t0 = time.time()
    YA = symmetric_vecs(m, at0=(uA0 + q) if wA else None)
    YB = free_vecs(m, at0=(uB0 + q) if wB else None)
    YD = free_vecs(m, at0=(uD0 + q) if wD else None)
    if not wB and not wD: YD = PC.canonical_translates(YD.astype(np.int8)).astype(np.int64)
    P = [pafs0(YA, h), 2 * pafs0(YB, h), pafs0(YD, h), np.zeros((1, h), np.int64)]
    tot1, sols = PC.join(P, h, K=a.maxy)
    log(f"v={v} m={m} q={q} layout {name} (u at w-positions: A {uA0}, B {uB0}, D {uD0}): stage 1 lists "
        f"{[len(YA), len(YB), len(YD)]}, triples {tot1} [{time.time() - t0:.0f}s]")
    if not tot1: return 0
    ZA = symmetric_vecs(m, at0=0) if wA else symmetric_vecs(m)
    ZB = free_vecs(m, zero_at=(0,)) if wB else free_vecs(m)
    ZD = free_vecs(m, zero_at=(0,)) if wD else free_vecs(m)
    found = 0
    for i, s in enumerate(sols):
        y = [YA[s[0]], YB[s[1]], YD[s[2]]]
        Wv = [np.zeros(m, np.int64) for _ in range(3)]
        for k, wk in enumerate((wA, wB, wD)):
            if wk: Wv[k][0] = 1
        u = [y[k] - q * Wv[k] for k in range(3)]
        c = (1, 2, 1)
        R = sum(c[k] * Wv[k][0] * np.array([u[k][j] + u[k][(-j) % m] for j in range(1, h + 1)]) for k in range(3))
        V = []
        for k, Z in enumerate((ZA, ZB, ZD)):
            V.append(c[k] * np.concatenate([pafs0(Z, h), cross(Z, u[k], h)], axis=1))
        V.append(np.zeros((1, 2 * h + 1), np.int64))
        target = np.concatenate([R, np.zeros(h + 1, np.int64)])
        tz, sz = PC.join(V, 2 * h + 1, K=a.nsol, target=target)
        for t in sz:
            Z = [ZA[t[0]], ZB[t[1]], ZD[t[2]]]
            seqs, M = lift(m, q, field, u, Wv, Z)
            ok = check(m, q, field, seqs)
            log(f"  triple {i}: {tz} z-solutions; lifted: symmetric Hadamard {ok}")
            if ok:
                found += 1
                rec = dict(n=v, method=f'propus_paley {name}', q=q, m=m, freeD=True, rows=[int(x.sum()) for x in seqs],
                           seqs=[''.join('+' if x > 0 else '-' for x in sq) for sq in seqs])
                if field:
                    k_ = field[1]
                    rec.update(group=[m] + [field[0]] * k_, automorphisms=[
                        np.block([[np.eye(1, dtype=np.int64), np.zeros((1, k_), np.int64)],
                                  [np.zeros((k_, 1), np.int64), M]]).tolist()])
                else:
                    def order(g):
                        k, x = 1, g % q
                        while x != 1: x = x * g % q; k += 1
                        return k
                    rec['gens'] = [next(g for g in range(1, v) if g % m == 1 and math.gcd(g, v) == 1 and order(g) == (q - 1) // 2)]
                with open(W.hits('sym_hits.jsonl'), 'a') as f: f.write(json.dumps(rec) + '\n')
        if i < 5 or i % 100 == 0: log(f"  triple {i}: z-solutions {tz} [{time.time() - t0:.0f}s]")
        if found: break
    log(f"v={v} layout {name}: {min(len(sols), i + 1)} of {tot1} triples tried, {found} verified [{time.time() - t0:.0f}s]")
    return found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('v', type=int); ap.add_argument('--q', type=int, default=0); ap.add_argument('--field', default='')
    ap.add_argument('--maxy', type=int, default=2000); ap.add_argument('--nsol', type=int, default=2)
    a = ap.parse_args()
    log = lambda s: print(s, flush=True)
    field = tuple(int(x) for x in a.field.split(',')) if a.field else None
    qs = [a.q] if a.q else [q for q in range(5, a.v + 1, 4) if a.v % q == 0 and math.gcd(q, a.v // q) == 1
                            and all(q % r for r in range(2, int(q ** 0.5) + 1))]
    for q in qs:
        m = a.v // q
        for lay in layouts(m, q):
            if search(a.v, m, q, field, lay, a, log): log(f"HIT v={a.v} q={q} layout {lay[0]}"); return


if __name__ == '__main__':
    main()

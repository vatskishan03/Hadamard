"""Skew Goethals-Seidel search on Paley-product structures: v = m q with q = 3 (mod 4) prime, q not dividing m,
and H = { h : h = 1 (mod m), h a square (mod q) } (order (q - 1)/2, odd, -1 not in H).

H fixes Z_m x {0} pointwise and acts on the Z_q factor as the quadratic residues, so it has 3m orbits:
{(a, 0)}, {a} x QR, {a} x NQR. Every block is X = f0 (x) 1 + f1 (x) Q + f2 (x) N with f0, f1, f2 in {+-1}^(Z_m),
and with Q* = N, QN = ((q+1)/4)(Q + N) + ((q-1)/2) the condition sum_k X_k X_k* = 4v becomes two equations in
the group ring Z[Z_m]:
    sum_k f0 f0* + ((q-1)/2)(f1 f1* + f2 f2*)                                   = 4v
    sum_k f0 f2* + f1 f0* + ((q+1)/4)(f1 f1* + f2 f2* + f2 f1*) + ((q-3)/4) f1 f2* = 0
a problem over the short group Z_m with an 8-letter alphabet per position. On v = 253 = 11 * 23 the Gaussian
models (diagonal or full covariance) predict 10^-2 to 10^-6 solutions in small low-energy lists, but the lists
hold thousands (solutions also come in clumps: B, C, D can each be translated by Z_m x {0}).

For each row-sum pattern the script builds low-energy lists (skew_gen.c) and all (or a sample of) the skew A,
counts every solution in the lists exactly with join4x, and verifies and logs the first few.

  python3 search/paley_product.py [--n N] [--steps S] [--amax N] [--maxo O] [--patterns K] v [v ...]
      -> hits in $HADAMARD_WORK/skew_hits.jsonl, one line per (v, pattern) in $HADAMARD_WORK/paley_product.jsonl
"""
import argparse, itertools, json, math, os, struct, subprocess, sys, time
import numpy as np
from pathlib import Path as _P
sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
import spence_lib as L
import skew_ktree as K
import verify


def is_prime(n):
    return n > 1 and all(n % p for p in range(2, int(n ** 0.5) + 1))


def structures(v, maxo):
    """(m, q, generator of H) for each admissible factorisation v = m q."""
    out = []
    for q in range(3, v + 1):
        if v % q or not is_prime(q) or q % 4 != 3: continue
        m = v // q
        if m % q == 0 or m < 3 or 3 * m > maxo: continue
        r = next(g for g in range(2, q) if all(pow(g, (q - 1) // f, q) != 1 for f in range(2, q) if (q - 1) % f == 0 and is_prime(f)))
        s = r * r % q                                    # generates the squares mod q
        h = next(x for x in range(1, v) if x % m == 1 and x % q == s)     # CRT lift: 1 mod m, s mod q
        out.append((m, q, h))
    return out


def patterns(v):
    """Row sums (b, c, d), odd, 0 < b <= c <= d, with 1 + b^2 + c^2 + d^2 = 4v."""
    res = []
    for b in range(1, int((4 * v) ** 0.5) + 1, 2):
        for c in range(b, int((4 * v) ** 0.5) + 1, 2):
            r = 4 * v - 1 - b * b - c * c
            if r < c * c: break
            d = int(round(r ** 0.5))
            if d * d == r and d % 2 == 1: res.append((b, c, d))
    return res


def vals16(st, Y, pm):
    out = np.empty((len(Y), st.C), np.int16)
    ch = max(1000, int(1.5e8 // (st.C * st.O)))
    for s in range(0, len(Y), ch): out[s:s + ch] = np.rint((st.class_values(Y[s:s + ch].astype(np.float64)) - pm) / 4)
    return out


def search(v, h, pat, n, steps, amax, seed, nsol=4):
    """Low-energy lists for structure <h> (h an int or a list of generators) and pattern pat, all skew A (or a
    sample of amax), exact count of the solutions in the lists (join4x), and the verified first few."""
    gens = list(h) if isinstance(h, (list, tuple)) else [h]
    G = L.AbGroup([v]); st = L.Structure(G, [G.perm_from_matrix(np.array([[g]])) for g in gens]); pm = v % 4
    rng = np.random.default_rng(seed)
    rep = np.array([np.flatnonzero(st.orb == o)[0] for o in range(st.O)]); partner = st.orb[G.neg()[rep]]
    pairs = [(o, partner[o]) for o in range(1, st.O) if o < partner[o]]
    if 2 ** len(pairs) <= amax:
        s = 1 - 2 * ((np.arange(2 ** len(pairs))[:, None] >> np.arange(len(pairs))) & 1)
    else:
        s = np.unique(rng.choice(np.array([-1, 1], np.int8), size=(int(amax * 1.05), len(pairs))), axis=0)[:int(amax)]
    YA = np.zeros((len(s), st.O), np.int8); YA[:, 0] = 1
    for j, (o, q) in enumerate(pairs): YA[:, o] = s[:, j]; YA[:, q] = -s[:, j]
    lists = {r: K.gen_pool(st, r, n, steps, rng, v) for r in sorted(set(pat))}
    Ys = [YA] + [lists[r] for r in pat]
    Vs = [vals16(st, Y, pm) for Y in Ys]
    sz = [len(V) for V in Vs]
    order = np.argsort(sz); perm = [order[0], order[3], order[1], order[2]]
    path = os.path.join(W.tmpdir(), f'pp_{v}_{os.getpid()}.bin')
    with open(path, 'wb') as f:
        f.write(struct.pack('iiiii', *[sz[k] for k in perm], st.C)); f.write(np.full(st.C, -pm, np.int16).tobytes())
        for k in perm: f.write(np.ascontiguousarray(Vs[k]).tobytes())
    r = subprocess.run([W.bin('join4x'), path, str(nsol)], capture_output=True, text=True, check=True)
    os.remove(path)
    total = int(r.stdout.split('TOTAL')[-1])
    found = []
    for line in r.stdout.splitlines():
        if not line.startswith('SOL'): continue
        ids = list(map(int, line.split()[1:]))
        q4 = [None] * 4
        for slot, k in enumerate(perm): q4[k] = ids[slot]
        seqs = [Ys[k][q4[k]].astype(np.int64)[st.orb] for k in range(4)]
        tot = verify.paf_sum(seqs)
        if not (tot[0] == 4 * v and not tot[1:].any()): continue
        H = verify.goethals_seidel(*map(verify.circulant, seqs)).astype(np.float64)
        had = np.array_equal(H @ H.T, 4 * v * np.eye(4 * v)); skew = np.array_equal(H + H.T, 2 * np.eye(4 * v))
        if had and skew:
            found.append(dict(n=v, gens=gens, pattern=[int(x.sum()) for x in seqs], hadamard=True, skew=True,
                              method='lowE_exact', seed=seed, n_lists=n, steps=steps,
                              seqs=[''.join('+' if x > 0 else '-' for x in sq) for sq in seqs]))
    return st, sz, total, found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('v', type=int, nargs='+')
    ap.add_argument('--n', type=float, default=1e5); ap.add_argument('--steps', type=int, default=3)
    ap.add_argument('--amax', type=float, default=2 ** 17); ap.add_argument('--maxo', type=int, default=60)
    ap.add_argument('--patterns', type=int, default=0, help='at most this many row-sum patterns per structure (0: all)')
    ap.add_argument('--seed', type=int, default=1); ap.add_argument('--stop', action='store_true', help='stop a v at its first hit')
    a = ap.parse_args()
    out = W.WORK / 'paley_product.jsonl'
    for v in a.v:
        for m, q, h in structures(v, a.maxo):
            pats = patterns(v)
            if a.patterns: pats = pats[:a.patterns]
            hit = False
            for pat in pats:
                t0 = time.time()
                try:
                    st, sz, total, found = search(v, h, pat, int(a.n), a.steps, a.amax, a.seed)
                except Exception as e:
                    print(f"v={v} m={m} q={q} pattern={pat}: error {e!r}", flush=True); continue
                rec = dict(v=v, m=m, q=q, gens=[h], O=st.O, C=st.C, pattern=list(pat), lists=sz, solutions=total,
                           verified=len(found), seed=a.seed, n=a.n, steps=a.steps, seconds=round(time.time() - t0))
                with open(out, 'a') as f: f.write(json.dumps(rec) + '\n')
                print(json.dumps(rec), flush=True)
                if found:
                    with open(W.hits(), 'a') as f:
                        for x in found: f.write(json.dumps(x) + '\n')
                    hit = True
                    if a.stop: break
            print(f"== v={v} m={m} q={q} H=<{h}>: {'HIT' if hit else 'no hit'}", flush=True)


if __name__ == '__main__':
    main()

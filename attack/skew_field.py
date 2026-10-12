"""Skew Goethals-Seidel families over the additive group of GF(p^k) = Z_p^k, invariant under an odd-order
subgroup of the multiplicative group (cyclotomy): H = <w^((q-1)/d)> for a primitive element w, d odd.

  python3 skew_field.py rank p k           # structures with their Gaussian estimate and exact work
  python3 skew_field.py run p k d b,c,d    # exact search (find4b) or k-tree when the pools are large
"""
import itertools, json, math, os, sys
import numpy as np
sys.path.insert(0, '/home/user/Hadamard/attack')
import spence_lib as L
import skew_rank as R
import skew_group as SG
HERE = '/home/user/Hadamard/attack'


def primitive_companion(p, k):
    """Companion matrix of a primitive polynomial of degree k over GF(p) (order p^k - 1)."""
    q = p ** k
    for coeffs in itertools.product(range(p), repeat=k):
        if coeffs[0] == 0: continue
        M = np.zeros((k, k), dtype=np.int64)
        M[1:, :-1] = np.eye(k - 1, dtype=np.int64); M[:, -1] = [(-c) % p for c in coeffs]
        X = M.copy(); order = 1; I = np.eye(k, dtype=np.int64)
        ok = True
        for f in [f for f in range(2, q) if (q - 1) % f == 0 and all(f % r for r in range(2, int(f ** 0.5) + 1))]:
            e = (q - 1) // f
            if np.array_equal(mpow(M, e, p), I): ok = False; break
        if ok and np.array_equal(mpow(M, q - 1, p), I): return M
    raise ValueError


def mpow(M, e, p):
    R_ = np.eye(len(M), dtype=np.int64); B = M % p
    while e:
        if e & 1: R_ = R_ @ B % p
        B = B @ B % p; e >>= 1
    return R_


def structure(p, k, d):
    q = p ** k; W = primitive_companion(p, k)
    g = mpow(W, (q - 1) // d, p)
    G = L.AbGroup([p] * k)
    return G, L.Structure(G, [G.perm_from_matrix(g)]), g


if __name__ == '__main__':
    p, k = int(sys.argv[2]), int(sys.argv[3]); q = p ** k
    if sys.argv[1] == 'rank':
        rng = np.random.default_rng(1)
        for d in [d for d in range(3, q, 2) if (q - 1) % d == 0]:
            G, st, g = structure(p, k, d)
            if not (8 <= st.O <= 120): continue
            res = R.rank_structure(st, q, 1.2e7, rng=rng)
            for lE2, lE, qq, lN in sorted(res, reverse=True)[:2]:
                NA = 10 ** lN[0]; Ns = sorted(10 ** x for x in lN[1:]); work = NA * Ns[2] + Ns[0] * Ns[1]
                print(f"{lE2:7.2f} GF({p}^{k}) d={d} O={st.O} C={st.C} rows=(1,{qq}) log10E={lE:.2f} "
                      f"log10N=({','.join(f'{x:.1f}' for x in lN)}) work={work:.1e}", flush=True)
            if not res: print(f"   none GF({p}^{k}) d={d} O={st.O}: no admissible pattern", flush=True)
    else:
        d = int(sys.argv[4]); pat = [int(x) for x in sys.argv[5].split(',')]
        G, st, g = structure(p, k, d)
        res, work = SG.exact(q, G, st, tuple(pat), cap=float(sys.argv[6]) if len(sys.argv) > 6 else 6e10)
        print(f"GF({p}^{k}) d={d} O={st.O} C={st.C} rows=(1,{pat}) work={work:.1e}: "
              f"{'skipped' if res is None else str(res[0]) + ' solutions'}", flush=True)
        if res:
            for had, skew, seqs in res[1]:
                print(f"  VERIFIED Hadamard={had} skew={skew}", flush=True)
                if had and skew:
                    with open(os.path.join(HERE, 'skew_group_hits.json'), 'a') as f:
                        f.write(json.dumps(dict(n=q, mods=[p] * k, structure=f"singer{d}", generator=g.tolist(), hadamard=True,
                                                skew=True, pattern=[1, *pat], seqs=[''.join('+' if x > 0 else '-' for x in s) for s in seqs])) + '\n')

"""Szekeres difference sets with multiplier symmetry: skew-Hadamard matrices of order 2(w + 1) from a group of order w.

Let G be an abelian group of odd order w and a, b two +-1 functions on G with
- a of skew type: a(0) = 1, a(-x) = -a(x) for x != 0;
- b symmetric: b(-x) = b(x), with b(0) = 1 and row sum 1;
- PAF_a(s) + PAF_b(s) = -2 for every s != 0 (a Legendre pair).
(In set language: a 2-{w; (w-1)/2, (w-1)/2; (w-3)/2} supplementary difference set with a skew and b symmetric,
the "Szekeres difference sets".) With A, B the group-developed matrices (A + A^T = 2I, B = B^T, A A^T + B B^T =
(2w+2)I - 2J, row sums 1) the bordered array

        [  1    1    e    e  ]
    H = [ -1    1    e   -e  ]
        [ -e^T -e^T  A    B  ]
        [ -e^T  e^T -B   A^T ]

is a skew-Hadamard matrix of order 2(w + 1): H H^T = (2w + 2) I and H + H^T = 2I.

Both blocks are taken constant on the orbits of a subgroup H of odd order of Aut(G) (multipliers for cyclic G,
GL(2, p) parts as in skew_group.py for split primes). Every sign choice is listed exactly on each side and the two
lists are joined on all autocorrelation classes, so a structure that yields nothing is ruled out.
Generic counting is against this family (about two free signs per class equation), so it is a search for
algebraic coincidences in symmetric structures, which is why it is exhaustive and cheap.

  python3 search/szekeres.py MAXPAIRS w [w ...]      -> hits in $HADAMARD_WORK/szekeres_hits.jsonl
"""
import json, sys, time
import numpy as np
from pathlib import Path as _P; sys.path[:0] = [str(_P(__file__).resolve().parent)]
import workdir as W
import spence_lib as L
from gs_sweep_subgroups import subgroups
import skew_group as SG


def bordered(G, a, b):
    w = G.v
    diff = G.index(G.coords[None, :, :] - G.coords[:, None, :])        # diff[g, h] = h - g
    A, B = a[diff], b[diff]
    e = np.ones((1, w), dtype=np.int64)
    top = np.block([[np.array([[1, 1], [-1, 1]]), np.vstack([np.hstack([e, e]), np.hstack([e, -e])])]])
    left = np.hstack([np.vstack([-e.T, -e.T]), np.vstack([-e.T, e.T])])
    body = np.block([[A, B], [-B, A.T]])
    return np.vstack([top, np.hstack([left, body])])


def rowhash(V, R):
    return (V.astype(np.uint64) * R[None, :]).sum(1)


def search(G, st, maxpairs, label):
    w = G.v
    rep = np.array([np.flatnonzero(st.orb == o)[0] for o in range(st.O)])
    partner = st.orb[G.neg()[rep]]
    pairs = [(o, int(partner[o])) for o in range(1, st.O) if o < partner[o]]
    n = len(pairs)
    if n > maxpairs or n == 0: return None
    bits = 1 - 2 * ((np.arange(2 ** n, dtype=np.int64)[:, None] >> np.arange(n)) & 1)
    YA = np.zeros((2 ** n, st.O)); YA[:, 0] = 1
    YB = np.zeros((2 ** n, st.O)); YB[:, 0] = 1
    sz = np.array([st.sizes[o] for o, _ in pairs])
    for j, (o, q) in enumerate(pairs):
        YA[:, o] = bits[:, j]; YA[:, q] = -bits[:, j]
        YB[:, o] = bits[:, j]; YB[:, q] = bits[:, j]
    YB = YB[(bits * sz).sum(1) == 0]                   # row sum of b: 1 + 2 sum size*t = 1
    if not len(YB): return dict(status='no b with row sum 1', pairs=n)
    def vals(Y):
        out = np.empty((len(Y), st.C), np.int64)
        step = max(1000, int(2e8 // max(1, st.C * st.O * st.O)))
        for s in range(0, len(Y), step): out[s:s + step] = np.rint(st.class_values(Y[s:s + step]))
        return out
    VA, VB = vals(YA), vals(YB)
    R = np.random.default_rng(7).integers(1, 2 ** 62, st.C, dtype=np.int64).astype(np.uint64) | np.uint64(1)
    hA = rowhash(-2 - VA, R); hB = rowhash(VB, R)
    common = np.intersect1d(hA, hB)
    sols = []
    if len(common):
        ia = np.flatnonzero(np.isin(hA, common)); ib = np.flatnonzero(np.isin(hB, common))
        byh = {}
        for j in ib: byh.setdefault(int(hB[j]), []).append(j)
        for i in ia:
            for j in byh.get(int(hA[i]), []):
                if np.array_equal(VA[i] + VB[j], np.full(st.C, -2)): sols.append((i, j))
    found = []
    for i, j in sols[:2]:
        a = YA[i].astype(np.int64)[st.orb]; b = YB[j].astype(np.int64)[st.orb]
        Hm = bordered(G, a, b).astype(np.float64)
        n2 = 2 * w + 2
        ok = np.array_equal(Hm @ Hm.T, n2 * np.eye(n2)) and np.array_equal(Hm + Hm.T, 2 * np.eye(n2))
        found.append(dict(ok=bool(ok), a=''.join('+' if x > 0 else '-' for x in a), b=''.join('+' if x > 0 else '-' for x in b)))
    return dict(status='found' if any(f['ok'] for f in found) else 'ruled out', pairs=n, nA=len(YA), nB=len(YB),
                solutions=len(sols), found=found)


def structures(w, mino=3):
    """Cyclic G = Z_w with odd-order multiplier subgroups; for split primes also the GL(2, p) structures."""
    G = L.AbGroup([w])
    for Hs, gens in subgroups(w).items():
        if len(Hs) % 2 == 0 or len(Hs) == 1: continue
        st = L.Structure(G, [G.perm_from_matrix(np.array([[g]])) for g in gens])
        yield f"Z{w}<{','.join(map(str, gens))}>|{len(Hs)}", G, st
    for name, G2, st2 in SG.structures(w, mino, 10 ** 6) or []:
        yield f"nc:{name}", G2, st2


def main():
    maxpairs = int(sys.argv[1])
    out = W.WORK / 'szekeres.jsonl'
    for w in map(int, sys.argv[2:]):
        t0 = time.time(); nst = 0; hit = False
        for label, G, st in structures(w):
            res = search(G, st, maxpairs, label)
            if res is None: continue
            nst += 1
            rec = dict(w=w, order=2 * (w + 1), structure=label, O=st.O, C=st.C, **{k: v for k, v in res.items() if k != 'found'})
            with open(out, 'a') as f: f.write(json.dumps(rec) + '\n')
            if res.get('status') == 'found':
                hit = True
                with open(W.hits('szekeres_hits.jsonl'), 'a') as f:
                    for fd in res['found']:
                        if fd['ok']: f.write(json.dumps(dict(w=w, order=2 * (w + 1), structure=label, a=fd['a'], b=fd['b'])) + '\n')
                print(json.dumps(rec), flush=True)
                break
        print(f"== w={w} (order {2 * (w + 1)}): {'HIT' if hit else 'no hit'} in {nst} structures [{time.time()-t0:.0f}s]", flush=True)


if __name__ == '__main__':
    main()

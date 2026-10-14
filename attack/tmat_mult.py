#!/usr/bin/env python3
"""T-matrices of order t from a multiplier structure, through the four-list split (attack/NOTES.md section 6).

T-matrices: circulant 0/+-1 matrices T_1..T_4 of order t with partitioning supports and sum PAF = t delta. All four
are constant on the orbits of a multiplier group H <= Z_t^*. Fix E = the orbits in supp T_1 u supp T_2; then
X_1 = T_1 + T_2 and X_3 = T_1 - T_2 are +-1 on E, X_2 = T_3 + T_4 and X_4 = T_3 - T_4 are +-1 off E, and
sum PAF(X_i) = 2 sum PAF(T_i). So for a fixed E the four lists are independent: each is every sign vector on its
orbits, and they are joined exactly on the classes of <H, -1> (numpy sort join). E is drawn at random until a
quadruple is found.

    python3 attack/tmat_mult.py t g [--tries N] [--seed S] [--out FILE]
"""
import argparse, itertools, sys
import numpy as np


def closure(gens, m):
    S, todo = {1}, [1]
    while todo:
        x = todo.pop()
        for g in gens:
            y = x * g % m
            if y not in S:
                S.add(y); todo.append(y)
    return sorted(S)


def paf_rows(X):
    F = np.fft.fft(X.astype(np.float64), axis=1)
    return np.rint(np.fft.ifft(F * np.conj(F), axis=1).real).astype(np.int64)


def lists_for(orbits, sel, t, reps):
    """All +-1 assignments on the selected orbits (others 0), as sequences and class values."""
    k = len(sel)
    signs = np.array(list(itertools.product((1, -1), repeat=k)), dtype=np.int64) if k else np.zeros((1, 0), np.int64)
    X = np.zeros((len(signs), t), np.int64)
    for j, o in enumerate(sel):
        X[:, orbits[o]] = signs[:, j:j + 1]
    return X, paf_rows(X)[:, reps]


def join(L1, L2, L3, L4):
    """(i, j, k, l) with L1[i] + L2[j] + L3[k] + L4[l] = 0 on every class: sort-intersection of the pair sums."""
    C = L1.shape[1]
    L1, L2, L3, L4 = (x.astype(np.int16) for x in (L1, L2, L3, L4))
    left = (L1[:, None, :] + L2[None, :, :]).reshape(-1, C)
    right = (-(L3[:, None, :] + L4[None, :, :])).reshape(-1, C)
    kl = np.ascontiguousarray(left).view(np.dtype((np.void, 2 * C))).ravel()
    kr = np.ascontiguousarray(right).view(np.dtype((np.void, 2 * C))).ravel()
    _, il, ir = np.intersect1d(kl, kr, return_indices=True)
    if len(il) == 0:
        return None
    i, j = divmod(int(il[0]), len(L2)); k, l = divmod(int(ir[0]), len(L4))
    return i, j, k, l


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("t", type=int); ap.add_argument("g", type=int, nargs="+")
    ap.add_argument("--tries", type=int, default=2000); ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out")
    a = ap.parse_args()
    t = a.t
    H = closure(a.g, t)
    L = sorted(set(H) | {(-h) % t for h in H})
    seen, orbits = set(), [[0]]
    for x in range(1, t):
        if x not in seen:
            o = sorted({x * h % t for h in H}); seen.update(o); orbits.append(o)
    seen, reps = set(), []
    for x in range(1, t):
        if x not in seen:
            o = {x * h % t for h in L}; seen.update(o); reps.append(x)
    print(f"t = {t}, |H| = {len(H)}: {len(orbits)} orbits, {len(reps)} classes", flush=True)
    rng = np.random.default_rng(a.seed)
    for attempt in range(a.tries):
        inE = np.zeros(len(orbits), bool)                     # a balanced split of the orbits
        inE[rng.permutation(len(orbits))[:len(orbits) // 2]] = True
        E = [o for o in range(len(orbits)) if inE[o]]; F = [o for o in range(len(orbits)) if not inE[o]]
        X1, V1 = lists_for(orbits, E, t, reps)
        X2, V2 = lists_for(orbits, F, t, reps)
        hit = join(V1, V1, V2, V2)
        if hit:
            i, j, k, l = hit
            x1, x3, x2, x4 = X1[i], X1[j], X2[k], X2[l]
            T = [(x1 + x3) // 2, (x1 - x3) // 2, (x2 + x4) // 2, (x2 - x4) // 2]
            P = sum(paf_rows(np.array([x]))[0] for x in T)
            ok = (np.count_nonzero(np.array(T), axis=0) == 1).all() and P[0] == t and not P[1:].any()
            rows = ["".join("+" if v > 0 else "-" if v < 0 else "0" for v in x) for x in T]
            print(f"attempt {attempt}: T-matrices of order {t} found, verified = {ok}, row sums {[int(x.sum()) for x in T]}")
            for r in rows:
                print(r)
            if ok and a.out:
                open(a.out, "w").write("\n".join(rows) + "\n")
            return
    print("no T-matrices found")


if __name__ == "__main__":
    main()

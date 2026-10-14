#!/usr/bin/env python3
"""Sampled Szöllősi (ssxx) search for structures whose lists are too large to enumerate (attack/ssxx.py).

For each row pattern the lists of A and B (symmetric, <H, -1>-orbits) are complete when they are small, and C, D
(H-orbits) are uniform samples of the given size (DP over the orbits, as in bnw.Struct.sample). The four lists are
joined exactly (attack/join4x.c). With pools P_k and samples S_k the expected number of hits is E * prod(S_k / P_k),
E the Gaussian estimate of the structure (reliable to one or two orders for prime N, attack/NOTES.md section 15);
the script prints both. Repeat with other seeds to sample again.

    python3 attack/ssxx_sample.py N g [g2] [--sub 1e6] [--seed S] [--rounds R] [--log FILE]
    python3 attack/ssxx_sample.py sweep N1 [N2 ...] [--sub 1e6] [--minexp 0] [--log FILE]
"""
import argparse, json, math, sys, time
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ssxx
from bnw import Struct, join4


def fast_sample(sz, r, K, rng):
    """K uniform sign vectors y on orbits of sizes sz with sum sz_i y_i = r: backward DP of path counts (float, fine
    for sampling), then one vectorised step per orbit for all K samples."""
    O = len(sz); tot = int(sum(sz)); off = tot
    tab = np.zeros((O + 1, 2 * tot + 1))
    tab[O, off] = 1.0
    for i in range(O - 1, -1, -1):
        s = int(sz[i]); nxt = tab[i + 1]
        cur = np.zeros_like(nxt)
        cur[s:] += nxt[:-s] if s else nxt      # y_i = +1: remaining sum shifts down by s
        cur[:-s if s else None] += nxt[s:]     # y_i = -1
        tab[i] = cur
    if tab[0, off + r] == 0:
        return np.zeros((0, O), np.int8)
    Y = np.empty((K, O), np.int8)
    rem = np.full(K, r, np.int64)
    for i in range(O):
        s = int(sz[i])
        cp = tab[i + 1][np.clip(rem - s + off, 0, 2 * tot)] * ((rem - s + off >= 0) & (rem - s + off <= 2 * tot))
        cm = tab[i + 1][np.clip(rem + s + off, 0, 2 * tot)] * ((rem + s + off >= 0) & (rem + s + off <= 2 * tot))
        plus = rng.random(K) * (cp + cm) < cp
        Y[:, i] = np.where(plus, 1, -1)
        rem = rem - np.where(plus, s, -s)
    return Y


def sample_list(st, kind, r, n, rng, batch=200000):
    """n uniform sign vectors on the orbits of `kind` with row sum r, and their class values."""
    Ys, Vs = [], []
    left = n
    sz = st.sizes(kind)
    while left > 0:
        k = min(batch, left)
        Y = fast_sample(sz, r, k, rng)
        Ys.append(Y); Vs.append(np.rint(st.values(kind, Y.astype(np.float64))).astype(np.int16))
        left -= k
    return np.concatenate(Ys), np.concatenate(Vs)


def run(N, gens, sub=1e6, seed=1, rounds=1, log=None, minest=-3.0, K=4):
    st = Struct(N, gens, plain=True)
    pats = [(ssxx.gauss(st, p), p) for p in ssxx.patterns(st)]
    pats = sorted([x for x in pats if x[0] > minest], key=lambda x: -x[0])
    print(f"N = {N}, H = <{','.join(map(str, gens))}> of order {len(st.H)}: O_A = {len(st.sizesA)}, "
          f"O_L = {len(st.sizesL)}, classes = {st.C}; {len(pats)} patterns above 10^{minest}", flush=True)
    rng = np.random.default_rng(seed)
    found = []
    for rd in range(rounds):
        for est, pat in pats:
            t0 = time.time()
            pools = [st.counts(k).get(r, 0) for k, r in zip(ssxx.KINDS, pat)]
            lists, vecs, frac = [], [], 0.0
            for kind, r, P in zip(ssxx.KINDS, pat, pools):
                if P <= sub:
                    V, M = ssxx.enum(st, kind, r)
                    lists.append(V); vecs.append(("mask", M))
                else:
                    Y, V = sample_list(st, kind, r, int(sub), rng)
                    lists.append(V); vecs.append(("vec", Y))
                    frac += math.log10(sub / P)
            exp = est + frac
            sizes = [len(x) for x in lists]
            order = ssxx.best_pairing(sizes)
            sols, total = join4([lists[k] for k in order], st.C, K, None, target=0)
            print(f"  round {rd}, rows {tuple(int(x) for x in pat)}: pools {[f'{p:.2g}' for p in pools]}, lists {sizes}, "
                  f"expected 10^{exp:.2f}, found {total} ({time.time() - t0:.0f}s)", flush=True)
            for sol in sols:
                idx = [0] * 4
                for pos, k in enumerate(order):
                    idx[k] = sol[pos]
                seqs = []
                for k, (kind, i) in enumerate(zip(ssxx.KINDS, idx)):
                    typ, data = vecs[k]
                    if typ == "mask":
                        y = np.array([-1 if (int(data[i]) >> j) & 1 else 1 for j in range(len(st.sizes(kind)))])
                    else:
                        y = data[i].astype(np.int64)
                    seqs.append(st.expand(kind, y))
                ssxx.complex_hadamard(seqs)
                ssxx.goethals_seidel(seqs)
                rec = {"N": N, "gens": list(gens), "subgroup_order": len(st.H), "rowsums": [int(x.sum()) for x in seqs],
                       "sequences": ["".join("+" if v > 0 else "-" for v in x) for x in seqs]}
                found.append(rec)
                if log:
                    with open(log, "a") as f:
                        f.write(json.dumps(rec) + "\n")
            if found:
                print(f"  verified: complex Hadamard matrix of order {2 * N}", flush=True)
                return found
    return found


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "sweep":
        ap = argparse.ArgumentParser()
        ap.add_argument("cmd"); ap.add_argument("N", type=int, nargs="+"); ap.add_argument("--sub", type=float, default=1e6)
        ap.add_argument("--log"); ap.add_argument("--minexp", type=float, default=0.0); ap.add_argument("--seed", type=int, default=1)
        a = ap.parse_args()
        sweep(a.N, a.sub, a.log, a.minexp, seed=a.seed)
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("N", type=int); ap.add_argument("g", type=int, nargs="+")
    ap.add_argument("--sub", type=float, default=1e6); ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--rounds", type=int, default=1); ap.add_argument("--log")
    ap.add_argument("--minest", type=float, default=-3.0)
    a = ap.parse_args()
    f = run(a.N, a.g, a.sub, a.seed, a.rounds, a.log, a.minest)
    print(f"N = {a.N}: {'FOUND' if f else 'nothing'}")



def fast_gauss(st, pat, K=2000, seed=1):
    """ssxx.gauss with the vectorised sampler (same model: pools x class densities, one class dropped)."""
    rng = np.random.default_rng(seed)
    lN, mus, vars_, lats = 0.0, [], [], []
    for kind, r in zip(ssxx.KINDS, pat):
        n = st.counts(kind).get(r, 0)
        if n == 0:
            return -math.inf
        lN += math.log10(n)
        Y = fast_sample(st.sizes(kind), r, K, rng).astype(np.float64)
        V = st.values(kind, Y)
        mus.append(V.mean(0)); vars_.append(V.var(0))
        Vi = np.rint(V).astype(np.int64)
        lat = []
        for c in range(st.C):
            u = np.unique(Vi[:, c])
            g = int(np.gcd.reduce(np.diff(u))) if len(u) > 1 else 0
            lat.append((g, int(u[0] % g) if g else int(u[0])))
        lats.append(lat)
    lp, dropped = 0.0, False
    for c in range(st.C):
        mu = sum(x[c] for x in mus); var = sum(x[c] for x in vars_)
        g, res = 0, 0
        for lat in lats:
            g = math.gcd(g, lat[c][0]); res += lat[c][1]
        if var < 1e-9:
            if abs(mu) > 1e-6:
                return -math.inf
            continue
        if g and (-res) % g:
            return -math.inf
        if not dropped:
            dropped = True; continue
        lp += -0.5 * math.log10(2 * math.pi * var) - mu ** 2 / (2 * var) / math.log(10) + math.log10(max(g, 1))
    return lN + lp


def plan_sampled(N, sub=1e6, maxo=62, maxpair=3e10):
    """Every (H, pattern) with -1 not in H, ranked by the expected number of hits inside lists capped at `sub`."""
    from bnw import subgroups
    rows = []
    for H, gens in subgroups(N).items():
        if (N - 1) in H:
            continue
        st = Struct(N, gens, plain=True)
        if len(st.sizesA) > maxo:
            continue
        for pat in ssxx.patterns(st):
            pools = [st.counts(k).get(r, 0) for k, r in zip(ssxx.KINDS, pat)]
            if min(pools) == 0 or any(p > sub and len(st.sizes(k)) > 62 for p, k in zip(pools, ssxx.KINDS)):
                continue
            if any(p > sub and False for p in pools):
                continue
            # complete lists are enumerated with masks (<= 62 orbits); sampled lists need no masks
            if any(p <= sub and len(st.sizes(k)) > 62 for p, k in zip(pools, ssxx.KINDS)):
                continue
            e = fast_gauss(st, pat)
            if e == -math.inf:
                continue
            sizes = [min(p, sub) for p in pools]
            exp = e + sum(math.log10(min(1.0, sub / p)) for p in pools)
            cost = ssxx.pair_cost(sizes, ssxx.best_pairing(sizes))
            if cost <= maxpair:
                rows.append((exp, cost, gens, len(H), pat))
    rows.sort(key=lambda r: -r[0])
    return rows


def sweep(Ns, sub=1e6, log=None, minexp=0.0, maxtries=3, seed=1):
    for N in Ns:
        rows = plan_sampled(N, sub)
        if not rows or rows[0][0] < minexp:
            print(f"N = {N}: no sampled structure with expected >= 10^{minexp} "
                  f"(best {rows[0][0]:.2f})" if rows else f"N = {N}: no structure", flush=True)
            continue
        done = False
        for exp, cost, gens, h, pat in rows[:maxtries]:
            if exp < minexp:
                break
            # adaptive sample size: the expected count grows like S^k (k sampled lists), the join cost like S
            st = Struct(N, gens, plain=True)
            pools = [st.counts(kk).get(r, 0) for kk, r in zip(ssxx.KINDS, pat)]
            k = sum(1 for p in pools if p > sub)
            S = sub
            if k and exp > math.log10(8):
                S = max(5e4, sub * 10 ** ((math.log10(8) - exp) / k))
                exp = exp + sum(math.log10(min(1.0, S / p)) - math.log10(min(1.0, sub / p)) for p in pools)
            print(f"N = {N}: H = <{','.join(map(str, gens))}> (order {h}), rows {tuple(int(x) for x in pat)}, "
                  f"expected 10^{exp:.2f} in samples of {S:.2g}", flush=True)
            f = run_one(N, gens, pat, S, seed, log)
            if f:
                print(f"N = {N}: FOUND", flush=True); done = True
                break
        if not done:
            print(f"N = {N}: nothing", flush=True)


def run_one(N, gens, pat, sub, seed, log, K=2):
    st = Struct(N, gens, plain=True)
    rng = np.random.default_rng(seed)
    pools = [st.counts(k).get(r, 0) for k, r in zip(ssxx.KINDS, pat)]
    lists, vecs = [], []
    for kind, r, P in zip(ssxx.KINDS, pat, pools):
        if P <= sub:
            V, M = ssxx.enum(st, kind, r)
            lists.append(V); vecs.append(("mask", M))
        else:
            Y, V = sample_list(st, kind, r, int(sub), rng)
            lists.append(V); vecs.append(("vec", Y))
    sizes = [len(x) for x in lists]
    order = ssxx.best_pairing(sizes)
    t0 = time.time()
    sols, total = join4([lists[k] for k in order], st.C, K, None, target=0)
    print(f"  lists {sizes}: {total} solutions ({time.time() - t0:.0f}s)", flush=True)
    found = []
    for sol in sols:
        idx = [0] * 4
        for pos, k in enumerate(order):
            idx[k] = sol[pos]
        seqs = []
        for k, (kind, i) in enumerate(zip(ssxx.KINDS, idx)):
            typ, data = vecs[k]
            y = np.array([-1 if (int(data[i]) >> j) & 1 else 1 for j in range(len(st.sizes(kind)))]) \
                if typ == "mask" else data[i].astype(np.int64)
            seqs.append(st.expand(kind, y))
        ssxx.complex_hadamard(seqs)
        ssxx.goethals_seidel(seqs)
        rec = {"N": N, "gens": list(gens), "subgroup_order": len(st.H), "rowsums": [int(x.sum()) for x in seqs],
               "sequences": ["".join("+" if v > 0 else "-" for v in x) for x in seqs]}
        found.append(rec)
        if log:
            with open(log, "a") as f:
                f.write(json.dumps(rec) + "\n")
    return found


if __name__ == "__main__":
    main()

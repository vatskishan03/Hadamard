"""shapes4b.py -- extension of the shape enumeration: PRODUCT multiplier groups H_K = H_1 x H_2 x ... (one cyclic unit subgroup per prime-power coordinate of K,
at least two non-trivial and not already cyclic as a tuple group), with sign homomorphism (twist) per generator.  Same exact filters as shapes4."""
import sys, json, math, itertools, time
import numpy as np
import os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'gs'))
from shapes4 import *


def coord_cyclic_subgroups(m):
    """cyclic subgroups <u> of Z_m^*: list of (u, order)"""
    units = [u for u in range(1, m) if math.gcd(u, m) == 1] if m > 1 else [1]
    seen = {}
    for u in units:
        S = frozenset(pow(u, k, m) for k in range(0, 2 * m))
        if S not in seen: seen[S] = u
    out = []
    for S, u in seen.items():
        o = 1; x = u % m
        while x != 1 % m: x = x * u % m; o += 1
        out.append((u, o))
    return out


def product_subgroups(Kmods, twist=True, maxsize=400):
    mods = list(Kmods)
    per = [coord_cyclic_subgroups(m) for m in mods]
    one = tuple([1] * len(mods))
    res = {}
    for combo in itertools.product(*per):
        nontriv = [i for i, (u, o) in enumerate(combo) if o > 1]
        if len(nontriv) < 2: continue
        orders = [combo[i][1] for i in nontriv]
        if all(math.gcd(a, b) == 1 for a, b in itertools.combinations(orders, 2)): continue      # cyclic: already enumerated
        size = int(np.prod(orders))
        if size > maxsize: continue
        gens = []
        for i in nontriv:
            g = [1] * len(mods); g[i] = combo[i][0]; gens.append(tuple(g))
        for signs in itertools.product((1, -1) if twist else (1,), repeat=len(gens)):
            if any(s == -1 and combo[i][1] % 2 for s, i in zip(signs, nontriv)): continue
            # closure with signs
            pairs = {(one, 1)}; frontier = [(one, 1)]
            while frontier:
                x, sx = frontier.pop()
                for g, sg in zip(gens, signs):
                    y = (tuple((a * b) % m for a, b, m in zip(x, g, mods)), sx * sg)
                    if y not in pairs: pairs.add(y); frontier.append(y)
            if len(pairs) != size: continue            # sign homomorphism inconsistent
            key = frozenset(pairs)
            if key not in res: res[key] = ([(g, sg) for g, sg in zip(gens, signs)], size)
    return list(res.values())


def enumerate_shapes_b(n, modes='abc', Orange=(4, 28), qkinds=('z3', 'z2z2'), twist=True, rmin=5):
    for mode in modes:
        if mode == 'a': tlist = [(0, None)]
        elif mode == 'b': tlist = [(1, 'b')]
        else: tlist = [(3, qk) for qk in qkinds]
        for (t, qk) in tlist:
            Gsz = n - t
            if Gsz < 6: continue
            if mode == 'c' and qk == 'z2z2' and n % 4 != 3: continue
            for r in odd_prime_powers_dividing(Gsz, rmin):
                N = Gsz // r
                if N < 4: continue
                for Km0 in abelian_groups_pp(N):
                    if len(Km0) < 2: continue
                    if mode == 'a': variants = [(Km0, None)]
                    elif mode == 'b': variants = [(Km0, ([], [()]))]
                    else: variants = abelian_orderings(Km0, qk)
                    for (Km, qs) in variants:
                        for hg, hs in product_subgroups(Km, twist=twist):
                            no, nf = orbit_count(Km, hg)
                            if no < Orange[0] or no > Orange[1]: continue
                            yield dict(mode=mode, n=n, r=r, K=Km, hg=[[list(u), ph] for u, ph in hg], qspec=qs, qkind=qk, t=t, orb=no, forced=nf, hsz=hs, N=N, product=True)


def main(out='candidates_b.jsonl'):
    ns = json.load(open('open_n.json'))
    f = open(out, 'w'); t0 = time.time(); cnt = 0
    for n in ns:
        for d in enumerate_shapes_b(n):
            hg = [(tuple(u), ph) for u, ph in d['hg']]
            try: sh = Shape4(d['K'], d['r'], hg, d['qspec'])
            except AssertionError: d['status'] = 'killed:phi-not-invariant'; f.write(json.dumps(d) + '\n'); continue
            d['Did'] = int(sh.Did); d['rowsq'] = int(sh.rowsq); d['nU'] = sh.nU; d['R'] = sh.R
            lo, hi = (-(-sh.Did // (sh.q + 1)), (sh.Did // (sh.q - 3) if sh.q > 3 else 4 * sh.N)); d['defwin'] = [lo, hi]
            if sh.Did < 0: d['status'] = 'killed:identity-negative'
            elif sh.kCD.fw > hi: d['status'] = 'killed:forced-defects-exceed'
            else:
                ok, why = dp_filter(sh); d['status'] = 'live' if ok else 'killed:' + why
            f.write(json.dumps(d, default=lambda o: o.item() if hasattr(o, 'item') else str(o)) + '\n'); cnt += 1
        f.flush()
    print('done', cnt, '%.0fs' % (time.time() - t0))

if __name__ == '__main__':
    main()

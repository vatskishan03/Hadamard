"""Biclique placements (X+, X-) on Z_p (|X+| = |X-| = d) whose Paley sign matrix [chi(y - x)] is switching-equivalent to J - 2P
(the sign class of BH(d,6) for d = 7; for d = 3 it is rank one, use --rank1). Uses the affine group x -> a x + b (a square) to fix
0, 1 in X+ ... no: enumerate all X+ containing 0 (translations), then all a in squares to dedupe at the end."""
import itertools, sys, pickle
import numpy as np
p = int(sys.argv[1]); d = (p - 1) // 4; rank1 = '--rank1' in sys.argv
qr = {x * x % p for x in range(1, p)}
chi = lambda x: 0 if x % p == 0 else (1 if x % p in qr else -1)
sq = sorted(qr)
found = set()
for rest in itertools.combinations(range(1, p), d - 1):
    Xp = (0,) + rest
    others = [y for y in range(p) if y not in Xp]
    vecs = {}
    for y in others:
        v = tuple(chi(y - x) for x in Xp)
        if v[0] < 0: v = tuple(-t for t in v)
        vecs.setdefault(v, []).append(y)
    if rank1:
        for v, ys in vecs.items():
            for Xm in itertools.combinations(ys, d):
                found.add((Xp, tuple(sorted(Xm))))
        continue
    # columns must be +-(a with coordinate i flipped), one column per i
    for a in itertools.product((1, -1), repeat=d):
        if a[0] < 0: continue
        cols = []
        for i in range(d):
            w = list(a); w[i] = -w[i]
            if w[0] < 0: w = [-t for t in w]
            ys = vecs.get(tuple(w))
            if not ys: break
            cols.append(ys)
        else:
            for Xm in itertools.product(*cols):
                if len(set(Xm)) == d: found.add((Xp, tuple(sorted(Xm))))
# close under translations (X+ contained 0 only) to get all placements
allp = set()
for Xp, Xm in found:
    for b in range(p):
        allp.add((tuple(sorted((x + b) % p for x in Xp)), tuple(sorted((y + b) % p for y in Xm))))
print(f"p={p} d={d}: biclique placements {len(allp)}")
pickle.dump(sorted(allp), open(f'kdd_{p}.pkl', 'wb'))

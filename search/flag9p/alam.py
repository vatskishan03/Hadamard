"""Alamouti-Paley ansatz for symmetric conference matrices of order 9p+1 (p = 1 mod 4 prime).
U = even residues 2,4,..,p-1; b-pairs {x, p+1-x}; S = one from each b-pair (d = (p-1)/4 elements);
S1 = {s in S: chi(s) chi(1-s) = +1}, S2 = rest. Need P1 = sum_{S1} chi(s) w^{c_s} t^s, P2 likewise,
PAF(P1) + PAF(P2) = d delta (periodic, Z_p), phases c in Z3."""
import itertools, sys
import numpy as np
p = int(sys.argv[1]); d = (p - 1) // 4
qr = {x * x % p for x in range(1, p)}
chi = lambda x: 0 if x % p == 0 else (1 if x % p in qr else -1)
bpairs = [(x, p + 1 - x) for x in range(2, (p + 1) // 2, 2)] + [(x, p + 1 - x) for x in range((p + 1) // 2 + 1, p, 2) if False]
# b-pairs within U: x even in [2, p-1], partner p+1-x (even); list each once
seen = set(); bpairs = []
for x in range(2, p, 2):
    y = p + 1 - x
    if x not in seen and y not in seen:
        bpairs.append((x, y)); seen |= {x, y}
assert len(bpairs) == d, (len(bpairs), d)
w = np.exp(2j * np.pi / 3)
nsol = 0
for pick in itertools.product((0, 1), repeat=d):
    S = [bp[k] for bp, k in zip(bpairs, pick)]
    S1 = [s for s in S if chi(s) * chi(1 - s) == 1]; S2 = [s for s in S if s not in S1]
    # necessary: difference multiplicities within S1 and within S2 never equal 1
    from collections import Counter
    cnt = Counter()
    for T in (S1, S2):
        for a in T:
            for b in T:
                if a != b: cnt[(a - b) % p] += 1
    if any(v == 1 for v in cnt.values()):
        continue
    # phases: fix first element of each part to 0
    free = [s for s in S1[1:]] + [s for s in S2[1:]]
    for ph in itertools.product(range(3), repeat=len(free)):
        c = {s: 0 for s in S}; c.update(dict(zip(free, ph)))
        ok = True
        for h in cnt:
            tot = 0
            for T in (S1, S2):
                for a in T:
                    b = (a - h) % p
                    if b in T and b != a:
                        tot += chi(a) * chi(b) * w ** ((c[a] - c[b]) % 3)
            if abs(tot) > 1e-9:
                ok = False; break
        if ok:
            nsol += 1
            if nsol <= 5:
                print("SOL S1", S1, "S2", S2, "phases", [c[s] for s in S1], [c[s] for s in S2])
print(f"p={p} d={d} solutions={nsol}")

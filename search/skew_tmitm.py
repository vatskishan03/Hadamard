"""Skew-compatible T-structures T = (I, X, -eps X^T, Y) over Z_t, invariant under a multiplier group H = <mu> (-1 not in H):
every orbit pair (O, -O) of nonzero elements carries X on O, X on -O (sign +-1), or Y = +-1 on O and -+1 on -O.
Condition 2 X X* + Y Y* = (t - 1) delta. Spectral meet in the middle over the choice of Y-pairs.
    python3 tmitm.py t mu [allowY=1]"""
import sys, itertools, numpy as np
t, mu = int(sys.argv[1]), int(sys.argv[2]); allowY = int(sys.argv[3]) if len(sys.argv) > 3 else 1
def horb(x):
    o = [x]; y = x * mu % t
    while y != x: o.append(y); y = y * mu % t
    return o
seen = set(); pairs = []
for x in range(1, t):
    if x in seen: continue
    o = horb(x); no = sorted({(-y) % t for y in o})
    if set(o) & set(no): sys.exit("H-orbit closed under -1: inadmissible")
    seen.update(o); seen.update(no); pairs.append((sorted(o), no))
P = len(pairs)
# characters: representatives of H-orbits on Z_t (dual identified with Z_t); include 0
cs = []; cseen = set()
for k in range(t):
    if k in cseen: continue
    cseen.update(horb(k) if k else [0]); cs.append(k)
W = np.exp(2j * np.pi * np.outer(np.arange(t), cs) / t)          # t x C
ind = lambda S: np.bincount(S, minlength=t).astype(float)
Opos = np.array([ind(o) @ W for o, no in pairs])                 # P x C : spectrum of indicator of O
Oneg = np.array([ind(no) @ W for o, no in pairs])
print(f"t = {t}, mu = {mu}, |H| = {len(horb(1))}, orbit pairs {P}, characters {len(cs)}", flush=True)
target = t - 1
hits = []
HV = np.random.default_rng(1).integers(1, 2**31, size=len(cs)).astype(np.int64)
for k in (range(P + 1) if allowY else [0]):
    for Ysel in itertools.combinations(range(P), k):
        rest = [i for i in range(P) if i not in Ysel]
        # Y spectra: sum_i s_i (Opos_i - Oneg_i), s in +-1, first sign fixed (global sign of Y is free)
        if k:
            YB = np.array([Opos[i] - Oneg[i] for i in Ysel])
            sg = np.array(list(itertools.product((1, -1), repeat=k - 1)))
            sg = np.hstack([np.ones((len(sg), 1)), sg]) if k > 1 else np.ones((1, 1))
            Yp = np.abs(sg @ YB) ** 2
        else:
            sg = np.zeros((1, 0)); Yp = np.zeros((1, len(cs)))
        need = (target - Yp) / 2.0                                  # required |X^|^2
        tab = {}
        NR = np.rint(need * 4).astype(np.int64)
        for r, v in enumerate(NR):
            tab.setdefault(v.tobytes(), []).append(r)
        tkeys = np.unique(NR @ HV)
        # X: each rest pair -> (side, sign): 4 options, first one's sign fixed
        m = len(rest)
        if m == 0:
            continue
        opts = np.array([Opos[i] for i in rest]), np.array([Oneg[i] for i in rest])
        # enumerate in chunks over 4^m (sign of first fixed -> 2*4^(m-1))
        total = 2 * 4 ** (m - 1)
        CH = 1 << 18
        for st in range(0, total, CH):
            idx = np.arange(st, min(total, st + CH))
            side = np.zeros((len(idx), m), int); sign = np.ones((len(idx), m))
            x = idx.copy()
            side[:, 0] = x % 2; x //= 2
            for j in range(1, m):
                c = x % 4; x //= 4
                side[:, j] = c % 2; sign[:, j] = np.where(c // 2, -1, 1)
            spec = np.zeros((len(idx), len(cs)), complex)
            for j in range(m):
                spec += sign[:, j:j + 1] * np.where(side[:, j:j + 1] == 0, opts[0][j], opts[1][j])
            pw = np.rint(np.abs(spec) ** 2 * 4).astype(np.int64)
            hk = pw @ HV
            cand = np.nonzero(np.isin(hk, tkeys))[0]
            for r in cand:
                lst = tab.get(pw[r].tobytes())
                if lst:
                    for yi in lst:
                        hits.append((Ysel, tuple(sg[yi]) if k else (), [rest[j] for j in range(m)], side[r].tolist(), sign[r].tolist()))
                        if len(hits) <= 3: print("  hit: Y pairs", Ysel, "signs", hits[-1][1], "X sides", side[r].tolist(), "signs", sign[r].tolist(), flush=True)
print(f"t = {t}, mu = {mu}: {len(hits)} solutions (up to the fixed signs)", flush=True)

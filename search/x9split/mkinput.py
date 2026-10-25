"""input files for sa_split: mode 'orb' (H-orbit slots, 19 shift reps) or 'free' (single positions, d = 1..95)"""
import sys, numpy as np
sys.path.insert(0, '.')
from dk764 import n, quad, slots, alpha
which, mode = int(sys.argv[1]), sys.argv[2]
Q = quad(which); sl = slots()
lines = [str(n)] + [' '.join(map(str, Q[k])) for k in range(4)]
if mode == 'orb':
    ids = [[k * 39 + int(sl[t]) for t in range(n)] for k in range(4)]
    fix = [(k * 39, k) for k in range(4)]
    seen = set(); dreps = []
    for d in range(1, n):
        if d in seen: continue
        dreps.append(d)
        for h in [1, 39, 184, 109, 49]: seen.add(d * h % n); seen.add((-d * h) % n)
else:
    ids = [[k * n + t for t in range(n)] for k in range(4)]
    fix = [(0, 0)]
    dreps = list(range(1, (n - 1) // 2 + 1))
lines += [' '.join(map(str, ids[k])) for k in range(4)]
lines += [str(len(fix))] + [f"{s} {l}" for s, l in fix]
lines += [str(len(dreps)), ' '.join(map(str, dreps))]
open(f"in_{which}_{mode}.txt", "w").write('\n'.join(lines) + '\n')

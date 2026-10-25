#!/usr/bin/env python3
"""genpool.py DIR TAB F0 NPOS THR TAU K TARGET [CHUNK_PROPOSALS=2e8] [SEED0=1] [MAXCHUNKS=100000]
Restart-safe flat-necklace pool generator: runs `flat gen2` chunk by chunk (seed = SEED0 + chunk index), each chunk writes
DIR/c<i>.bin (distinct canonical necklaces of that chunk, raw uint64) and DIR/c<i>.done; stops when the distinct total
(over all done chunks) >= TARGET.  Logs cumulative distinct and CPU seconds to DIR/log.txt.  Merge: mergepool.py."""
import sys, os, subprocess, time, numpy as np
d, tab, f0, npos, thr, tau, K, target = sys.argv[1:9]
chunk = float(sys.argv[9]) if len(sys.argv) > 9 else 2e8
seed0 = int(sys.argv[10]) if len(sys.argv) > 10 else 1
maxch = int(sys.argv[11]) if len(sys.argv) > 11 else 100000
target = int(float(target)); Kf = float(K)
os.makedirs(d, exist_ok=True)
here = os.path.dirname(os.path.abspath(__file__))
nrest = int(chunk / Kf)
allu = None
def load_done():
    arrs = []
    for fn in sorted(os.listdir(d)):
        if fn.endswith('.done'):
            arrs.append(np.fromfile(os.path.join(d, fn[:-5] + '.bin'), dtype=np.uint64))
    return np.unique(np.concatenate(arrs)) if arrs else np.zeros(0, np.uint64)
cur = load_done(); cpu = 0.0
logf = open(os.path.join(d, 'log.txt'), 'a')
i = 0
while len(cur) < target and i < maxch:
    if os.path.exists(os.path.join(d, f'c{i}.done')): i += 1; continue
    fn = os.path.join(d, f'c{i}.bin')
    if os.path.exists(fn): os.remove(fn)
    t0 = time.time()
    out = subprocess.run([os.path.join(here, 'flat2'), 'gen2', tab, f0, npos, thr, tau, K, str(nrest), str(seed0 + i), fn],
                         capture_output=True, text=True)
    if out.returncode: print(out.stderr); sys.exit(1)
    nprop, nd, dt = out.stdout.split(); cpu += float(dt)
    new = np.fromfile(fn, dtype=np.uint64)
    before = len(cur); cur = np.unique(np.concatenate([cur, new]))
    open(os.path.join(d, f'c{i}.done'), 'w').write(f'{nprop} {nd} {dt}\n')
    logf.write(f'chunk {i} cpu {dt} chunk_distinct {nd} cum_distinct {len(cur)} new {len(cur)-before}\n'); logf.flush()
    i += 1
cur.tofile(os.path.join(d, 'pool.bin'))
print('pool', len(cur), 'cpu this run %.1f' % cpu)

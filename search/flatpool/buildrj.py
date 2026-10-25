#!/usr/bin/env python3
"""buildrj.py OUT TAB C T spec_A spec_B spec_C spec_D [seed]   spec = n:pool.bin:f0:rot
Writes the rj.c input (int8 base keys) and OUT.map.npz (necklace masks per role, in list order)."""
import sys, subprocess, numpy as np, os
out, tab, C, T = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
specs = sys.argv[5:9]; seed = int(sys.argv[9]) if len(sys.argv) > 9 else 0
here = os.path.dirname(os.path.abspath(__file__)); rng = np.random.default_rng(seed)
keys = []; rots = []; maps = {}; ns = []
for role, sp in zip("ABCD", specs):
    n, path, f0, rot = sp.split(':'); n = int(float(n))
    masks = np.fromfile(path, dtype=np.uint64)
    sel = rng.permutation(len(masks))[:n] if n < len(masks) else np.arange(len(masks))
    m = masks[sel]; m.tofile(out + f'.{role}.masks')
    subprocess.run([os.path.join(here, 'flat2'), 'keys', tab, f0, out + f'.{role}.masks', out + f'.{role}.keys'], check=True, capture_output=True)
    K = np.fromfile(out + f'.{role}.keys', dtype=np.int16).reshape(-1, C)
    assert K.min() > -127 and K.max() < 127
    keys.append(K.astype(np.int8)); rots.append(int(rot)); maps[role] = m; ns.append(len(m))
    os.remove(out + f'.{role}.keys')
    print(role, 'necklaces', len(m), 'rot', rot, 'key range', K.min(), K.max(), flush=True)
with open(out, 'wb') as f:
    np.array([C, T] + rots + ns, dtype=np.int32).tofile(f)
    for K in keys: np.ascontiguousarray(K).tofile(f)
np.savez(out + '.map.npz', **maps)

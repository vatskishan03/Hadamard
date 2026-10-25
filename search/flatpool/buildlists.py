#!/usr/bin/env python3
"""buildlists.py OUT.bin TAB C nA:poolA nB:poolB nC:poolC nD:poolD [seed]
pool = path to pool.bin (canonical necklace masks) + .keys (made here by `flat keys`; f0 given as pool:f0).
Roles A: necklaces only (rotation 0).  B,C,D: every necklace x all C rotations (entry index i*C+k, key[(e+k)%C]).
Writes kt4 input and OUT.map.npz (masks of the used necklaces per role)."""
import sys, subprocess, numpy as np, os
out, tab, C = sys.argv[1], sys.argv[2], int(sys.argv[3])
specs = sys.argv[4:8]; seed = int(sys.argv[8]) if len(sys.argv) > 8 else 0
here = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(seed)
lists = []; maps = {}
for role, sp in zip("ABCD", specs):
    n, path, f0 = sp.split(':'); n = int(float(n))
    masks = np.fromfile(path, dtype=np.uint64)
    sel = rng.permutation(len(masks))[:n] if n < len(masks) else np.arange(len(masks))
    m = masks[sel]; m.tofile(out + f'.{role}.masks')
    subprocess.run([os.path.join(here, 'flat2'), 'keys', tab, f0, out + f'.{role}.masks', out + f'.{role}.keys'], check=True, capture_output=True)
    K = np.fromfile(out + f'.{role}.keys', dtype=np.int16).reshape(-1, C)
    maps[role] = m
    if role == 'A':
        V = K
    else:
        idx = (np.arange(C)[None, :] + np.arange(C)[:, None]) % C     # [k][e] -> (e+k)%C
        V = K[:, idx].reshape(-1, C)
    lists.append(V)
    print(role, 'necklaces', len(m), 'entries', len(V), flush=True)
T = np.full(C, -3, dtype=np.int16)
with open(out, 'wb') as f:
    np.array([len(v) for v in lists] + [C], dtype=np.int32).tofile(f)
    T.tofile(f)
    for v in lists: np.ascontiguousarray(v, dtype=np.int16).tofile(f)
np.savez(out + '.map.npz', **maps)

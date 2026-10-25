#!/usr/bin/env python3
"""decode.py v ulen MAP.npz F0A,F0B,F0C,F0D 'SOL x0 k0 x1 k1 x2 k2 x3 k3' [out.json] -> exact verification (PAF sum, GS array, H H^T = 4v I)"""
import sys, json, numpy as np
import gsflat, gsverify
v, ulen = int(sys.argv[1]), int(sys.argv[2])
mp = np.load(sys.argv[3]); f0s = [int(x) for x in sys.argv[4].split(',')]
tok = sys.argv[5].split()[1:]
S = gsflat.build(v, ulen); m = S['m']
def rotm(x, k): k %= m; return ((x << k) | (x >> (m - k))) & ((1 << m) - 1) if k else x
seqs = []; signs = []
for i, role in enumerate("ABCD"):
    x, k = int(tok[2 * i]), int(tok[2 * i + 1])
    mask = int(mp[role][x])
    mk = rotm(mask, -k)          # entry (x,k) has key[(e+k)%C] <-> mask rotated by -k
    f, s = gsverify.seq_from_mask(S, f0s[i], mk)
    seqs.append(f); signs.append(s)
res = gsverify.verify(seqs)
print(res)
if len(sys.argv) > 6 and all(bool(x) for k, x in res.items() if k != "v" and k != "order"):
    json.dump(dict(v=v, order=4 * v, U_order=ulen, rho=int(S['rho']), f0=f0s, orbit_signs_by_coset_index=[[int(a) for a in s] for s in signs],
                   row_sums=[int(f.sum()) for f in seqs], sequences=[''.join('+' if a > 0 else '-' for a in f) for f in seqs], verification={k: bool(x) for k, x in res.items()}), open(sys.argv[6], 'w'))

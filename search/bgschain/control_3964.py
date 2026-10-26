#!/usr/bin/env python3
"""Positive control: take an engine solution (G = Z2 x Z2 x Z13 x Z19 coordinates) for the 3964 pattern, convert it to the
repository certificate format (G = Z2 x Z2 x Z247, z = CRT(z13, z19)), re-verify the PAF conditions, build the bordered
Goethals-Seidel matrix with the border of hadamard/3964/certificate.json (verify.build_bordered_gs) and check H H^T = 3964 I exactly.
    python3 control_3964.py SOLUTIONS.jsonl [line]"""
import sys, os, json, copy
import numpy as np
ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))      # repository root
sys.path.insert(0, os.path.join(ROOT, 'verify'))
import verify as V

def main(sol_path, idx=0):
    cert = json.load(open(os.path.join(ROOT, 'hadamard', '3964', 'certificate.json')))
    recs = [json.loads(l) for l in open(sol_path)]
    r = recs[idx]
    assert r['group'] == [2, 2, 13, 19]
    mods4 = [2, 2, 13, 19]; mods3 = [2, 2, 247]
    v = 988
    c4 = np.array(np.unravel_index(np.arange(v), mods4)).T
    # CRT: z in Z247 with z = z13 mod 13, z = z19 mod 19
    z = np.array([next(zz for zz in range(247) if zz % 13 == a and zz % 19 == b) for a, b in zip(c4[:, 2], c4[:, 3])])
    idx3 = np.ravel_multi_index((c4[:, 0], c4[:, 1], z), mods3)
    seqs3 = []
    for s in r['sequences']:
        x4 = np.array([1 if ch == '+' else -1 for ch in s])
        x3 = np.zeros(v, dtype=np.int64); x3[idx3] = x4
        seqs3.append("".join('+' if t > 0 else '-' for t in x3))
    c = copy.deepcopy(cert); c['sequences'] = seqs3
    c['rowsums'] = [int(sum(1 if ch == '+' else -1 for ch in s)) for s in seqs3]
    assert c['rowsums'] == cert['rowsums'], (c['rowsums'], cert['rowsums'])
    print('same as certificate sequences:', [a == b for a, b in zip(seqs3, cert['sequences'])])
    H, seqs = V.build_bordered_gs(c)
    H = H.astype(np.float32)
    n = H.shape[0]
    G = H @ H.T
    ok = np.array_equal(G, 3964 * np.eye(n, dtype=np.float32)) and n == 3964 and set(np.unique(H)) <= {-1.0, 1.0}
    print('order', n, 'H H^T == 3964 I:', ok)
    return ok

if __name__ == '__main__':
    sys.exit(0 if main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 0) else 1)

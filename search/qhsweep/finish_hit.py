#!/usr/bin/env python3
"""From a hit (a line of hits_<mode>_<N>.jsonl, written by sched.py or sliced.py) to a certificate and the matrix of order
4(2N + 1): the certificate C<2N>.json in the format of complex-hadamard/ (construction szollosi-ssxx, or
szollosi-bordered for a hit with a border) is written to --out, rebuilt from the file alone by verify/verify_complex.py
and checked exactly (K K* = 2N I in integers, entries 1, -1, i, -i, SHA-256 of K); then verify.complex_miyamoto(q, Kr, Ki)
for q = 2N + 1 (a prime = 3 mod 4) gives the quaternary Hadamard matrix of order 4q, and H H* = 4q I is checked
blockwise (float32 products of 0, +-1 entries: every partial sum is an integer below 2^24, hence exact).
Nothing is written into the repository: copy the certificate to complex-hadamard/ by hand.

    python3 search/qhsweep/finish_hit.py HITS.jsonl [--line 0] [--out DIR] [--skip-miyamoto]
"""
import argparse, hashlib, json, math, sys, time
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import qhlib as Q
sys.path.insert(0, str(Q.ROOT / 'verify'))
import verify, verify_complex


def certificate(hit):
    N = hit['N']
    if 'border' in hit:
        return {'order': 2 * N, 'construction': 'szollosi-bordered', 'N': N, 'group': [N - 1], 'symmetry': hit['gens'],
                'rowsums': hit['rowsums'], 'border': hit['border'], 'sequences': hit['sequences']}
    return {'order': 2 * N, 'construction': 'szollosi-ssxx', 'N': N, 'multipliers': hit['gens'],
            'subgroup_order': hit['subgroup_order'], 'rowsums': hit['rowsums'], 'sequences': hit['sequences']}


def blockwise_check(q, Hr, Hi, bs=1024):
    """H H* = 4q I for H = Hr + i Hi of order 4q, in row blocks"""
    n = 4 * q
    R, I = Hr.astype(np.float32), Hi.astype(np.float32)
    ok = True
    for s in range(0, n, bs):
        re = R[s:s + bs] @ R.T + I[s:s + bs] @ I.T
        im = I[s:s + bs] @ R.T - R[s:s + bs] @ I.T
        re[np.arange(re.shape[0]), s + np.arange(re.shape[0])] -= n
        ok = ok and not re.any() and not im.any()
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('hits'); ap.add_argument('--line', type=int, default=0)
    ap.add_argument('--out', default=str(Q.QW / 'cert')); ap.add_argument('--skip-miyamoto', action='store_true')
    a = ap.parse_args()
    cert = certificate(json.loads(Path(a.hits).read_text().splitlines()[a.line]))
    Kr, Ki = verify_complex.build(cert)
    verify_complex.check(Kr, Ki)                                   # raises unless K K* = 2N I exactly, entries units
    cert['sha256'] = hashlib.sha256(verify_complex.text(Kr, Ki).encode()).hexdigest()
    out = Path(a.out) / f"C{cert['order']}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(cert))
    c2 = json.loads(out.read_text())                               # rebuilt from the certificate file alone
    Kr, Ki = verify_complex.build(c2)
    verify_complex.check(Kr, Ki)
    assert hashlib.sha256(verify_complex.text(Kr, Ki).encode()).hexdigest() == c2['sha256']
    print(f"certificate {out}: K K* = {c2['order']} I exact, sha256 {c2['sha256']}", flush=True)
    q = c2['order'] + 1
    if a.skip_miyamoto:
        return
    if q % 4 != 3 or any(q % p == 0 for p in range(2, math.isqrt(q) + 1)):
        print(f'q = {q} is not a prime = 3 (mod 4): no complex Miyamoto matrix'); return
    t = time.time()
    Hr, Hi = verify.complex_miyamoto(q, np.asarray(Kr, dtype=np.int64), np.asarray(Ki, dtype=np.int64))
    assert Hr.shape == (4 * q, 4 * q) and np.all(np.abs(Hr.astype(np.int16)) + np.abs(Hi.astype(np.int16)) == 1)
    ok = blockwise_check(q, Hr, Hi)
    print(f"complex_miyamoto(q = {q}): H of order {4 * q}, entries 1, -1, i, -i, H H* = {4 * q} I {'holds' if ok else 'FAILS'} ({time.time() - t:.0f}s)", flush=True)
    if not ok:
        sys.exit(1)


if __name__ == '__main__':
    main()

"""Turn 'HIT verified A <+-> B <+-> D <+->' lines of p65 / pz5 / prow into sym_hits.jsonl records.

    python3 attack/hit2rec.py LOG n gens [q]

n is the group order (Z_n), gens the multiplier generators (comma separated) and q the Z_q-border (0 = none).
The matrix is rebuilt (propus.symmetric_matrix_v3 for q = 0, propus.qbordered_matrix for q >= 1), checked
(H = H^T, H H^T = N I) and the record appended to $HADAMARD_WORK/sym_hits.jsonl; then run attack/sym_certify.py
(q = 0) or attack/border_certify.py (q >= 1).
"""
import json, sys
import numpy as np
from pathlib import Path as _P
sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import workdir as W
import propus as PR
import verify


def main():
    log, n, gens = sys.argv[1], int(sys.argv[2]), [int(x) for x in sys.argv[3].split(',')]
    q = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    seen = set()
    for line in open(log):
        if not line.startswith('HIT verified'): continue
        f = line.split()
        A, B, D = (np.array([1 if c == '+' else -1 for c in f[i]], dtype=np.int64) for i in (3, 5, 7))
        key = (f[3], f[5], f[7])
        if key in seen: continue
        seen.add(key)
        tot = verify.paf_sum([A, B, B, D])
        if q:
            tgt = np.zeros(n, np.int64); tgt[::q] = -4 * q
            assert tot[0] == 4 * n and np.array_equal(tot[1:], tgt[1:]), 'wrong PAF targets'
            Hm, R = PR.qbordered_matrix(A, B, D, q) if q > 1 else PR.bordered_matrix(A, B, D)
            N = 4 * (n + q)
        else:
            assert tot[0] == 4 * n and not tot[1:].any(), 'wrong PAF sum'
            Hm, R = PR.symmetric_matrix_v3(A, B, D), None
            N = 4 * n
        ok = Hm is not None and np.array_equal(Hm, Hm.T) and np.array_equal(Hm @ Hm.T, N * np.eye(N, dtype=Hm.dtype))
        print(f'order {N}: symmetric Hadamard {ok}')
        if not ok: continue
        rec = dict(n=n, gens=gens, symB=False, freeD=True, plain=False, rows=[int(A.sum()), int(B.sum()), int(D.sum())],
                   seqs=[''.join('+' if x > 0 else '-' for x in s) for s in (A, B, D)])
        if q == 1: rec.update(border=True, order=N, K=R.tolist())
        elif q: rec.update(border=True, q=q, order=N, R=R.tolist())
        with open(W.hits('sym_hits.jsonl'), 'a') as fh: fh.write(json.dumps(rec) + '\n')
        print('recorded', rec['rows'])


if __name__ == '__main__':
    main()

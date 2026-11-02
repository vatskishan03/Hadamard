"""Conversion of a hit (A, B, D on Z_v) into the symmetric Hadamard matrix of order 4v and a repository certificate.

The matrix is the second propus array of search/NOTES.md section 17 (the form of the repository's orders 532 and 876):

        [ A    B    B    D  ]
    H = [ B^T  DR  -A   -BR ]            A, B, D circulant with first rows the sequences, R the back identity,
        [ B^T -A  -DR    BR ]            XR symmetric for circulant X, A symmetric  =>  H = H^T.
        [ D^T -BR   BR  -A  ]

`build_H` below is written here (not imported); `write_certificate` produces certificate.json and sequences.txt in the
repository's format (same keys, same order, indent 1) with the SHA-256 of the canonical CSV of H.  The repository's own
verify/verify.py accepts the result (tested on 532 and 876, p3cert_test.py, and on the pipeline's own output).

  python3.11 p3cert.py --v 239 --gens 10 --hits hits.txt --out DIR [--matrix]
      hits.txt: lines "maskA maskB maskD" (hex orbit masks printed by p3 -hits)
"""
import argparse, gzip, hashlib, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from p3struct import Structure
import p3run as P


def circ(x):
    x = np.asarray(x, dtype=np.int64)
    n = len(x)
    i = np.arange(n)
    return x[(i[None, :] - i[:, None]) % n]          # C[i, j] = x[(j - i) mod n]


def build_H(A, B, D):
    n = len(A)
    cA, cB, cD = circ(A), circ(B), circ(D)
    R = np.eye(n, dtype=np.int64)[::-1]
    BR, DR = cB @ R, cD @ R
    return np.block([[cA, cB, cB, cD], [cB.T, DR, -cA, -BR], [cB.T, -cA, -DR, BR], [cD.T, -BR, BR, -cA]])


def canonical_csv_digest(H):
    """SHA-256 of numpy.savetxt(H, fmt='%d', delimiter=',') bytes (the repository's canonical encoding)."""
    h = hashlib.sha256()
    for i in range(0, H.shape[0], 256):
        blk = H[i:i + 256]
        lines = [(','.join('1' if x > 0 else '-1' for x in row) + '\n') for row in blk]
        h.update(''.join(lines).encode())
    return h.hexdigest()


def neg_reps(st, x):
    """One representative (the smallest element) of every H-orbit of Z_v on which the sequence is -1."""
    reps = []
    for i, r in enumerate(st.B.reps):
        if x[r] < 0:
            reps.append(int(r))
    return sorted(reps)


def make_certificate(st, gens, A, B, D):
    v = st.v
    A, B, D = (np.asarray(x, dtype=np.int64) for x in (A, B, D))
    assert P.exact_check(st, A, B, D), 'not a solution'
    H = build_H(A, B, D)
    n = 4 * v
    assert np.array_equal(H, H.T), 'H is not symmetric'
    G = H @ H.T
    assert np.array_equal(G, n * np.eye(n, dtype=np.int64)), 'H H^T != 4v I'
    cert = {
        'order': n, 'construction': 'propus', 'symmetric': True, 'group': f'Z_{v}', 'v': v,
        'multipliers': list(gens), 'subgroup_order': len(st.H),
        'rowsums': [int(A.sum()), int(B.sum()), int(D.sum())],
        'blocks': [neg_reps(st, x) for x in (A, B, D)],
        'sha256': canonical_csv_digest(H),
    }
    return cert, H


def dump_cert(cert):
    """certificate.json text in the repository's layout (indent 1, short lists on one line, one block per line)."""
    def lst(a):
        return '[' + ', '.join(str(x) for x in a) + ']'
    lines = ['{']
    keys = list(cert)
    for i, k in enumerate(keys):
        comma = ',' if i < len(keys) - 1 else ''
        val = cert[k]
        if k == 'blocks':
            lines.append(' "blocks": [')
            lines.extend('  ' + lst(b) + (',' if j < len(val) - 1 else '') for j, b in enumerate(val))
            lines.append(' ]' + comma)
        elif isinstance(val, list):
            lines.append(f' {json.dumps(k)}: {lst(val)}{comma}')
        else:
            lines.append(f' {json.dumps(k)}: {json.dumps(val)}{comma}')
    lines.append('}')
    return '\n'.join(lines) + '\n'


def seq_line(x):
    return ''.join('+' if t > 0 else '-' for t in x)


def write_certificate(outdir, st, gens, A, B, D, matrix=False):
    cert, H = make_certificate(st, gens, A, B, D)
    d = Path(outdir)
    d.mkdir(parents=True, exist_ok=True)
    (d / 'certificate.json').write_text(dump_cert(cert))
    (d / 'sequences.txt').write_text('\n'.join(seq_line(x) for x in (A, B, D)) + '\n')
    if matrix:
        with gzip.open(d / f'H{4 * st.v}.csv.gz', 'wb', compresslevel=9) as f:
            for i in range(0, H.shape[0], 256):
                f.write(''.join((','.join('1' if x > 0 else '-1' for x in row) + '\n') for row in H[i:i + 256]).encode())
    return cert


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--v', type=int, required=True)
    ap.add_argument('--gens', required=True, help='comma-separated multipliers generating H')
    ap.add_argument('--hits', required=True)
    ap.add_argument('--out', required=True, help='directory; certificates go to OUT/symmetric/<order>[_k]/')
    ap.add_argument('--matrix', action='store_true')
    a = ap.parse_args()
    gens = [int(x) for x in a.gens.split(',')]
    st = Structure(a.v, gens)
    seen = set()
    n = 0
    for line in Path(a.hits).read_text().split('\n'):
        if not line.strip():
            continue
        ma, mb, md = (int(x, 16) for x in line.split()[:3])
        if (ma, mb, md) in seen:
            continue
        seen.add((ma, mb, md))
        A, B, D = P.masks_to_seqs(st, ma, mb, md)
        d = Path(a.out) / 'symmetric' / (str(4 * a.v) if n == 0 else f'{4 * a.v}_{n}')
        cert = write_certificate(d, st, gens, A, B, D, a.matrix)
        print(f'wrote {d}  rows {cert["rowsums"]}  sha256 {cert["sha256"]}')
        n += 1
    print(f'{n} certificate(s)')


if __name__ == '__main__':
    main()

"""Turn a "SOL A ... B ..." line of spin_ls.c into a verified hit in $HADAMARD_WORK/skew_hits.jsonl.

  python3 attack/spin_sol.py v h mu SOLFILE

spin_ls can be run directly on an input file written by spin_skew.py (csearch) so that a long run is not tied
to a Python process; this script rebuilds the four blocks, checks H H^T = 4vI and H + H^T = 2I, and logs the hit.
"""
import json, sys
from pathlib import Path as _P
sys.path[:0] = [str(_P(__file__).resolve().parent), str(_P(__file__).resolve().parent.parent / 'verify')]
import numpy as np
import spin_skew as SS
import workdir as W
import verify


def main():
    v, h, mu = map(int, sys.argv[1:4])
    lines = [l for l in open(sys.argv[4]) if l.startswith('SOL')]
    if not lines: sys.exit("no SOL line")
    tok = lines[0].split(); ia = tok.index('A'); ib = tok.index('B')
    A = np.array(tok[ia + 1:ib], dtype=float); B = np.array(tok[ib + 1:], dtype=float)
    sp = SS.Spin(v, [h], mu)
    seqs = sp.sequences(A, B)
    tot = verify.paf_sum(seqs)
    ok = tot[0] == 4 * v and not tot[1:].any()
    H = verify.goethals_seidel(*map(verify.circulant, seqs)).astype(np.float64)
    had = np.array_equal(H @ H.T, 4 * v * np.eye(4 * v)); skew = np.array_equal(H + H.T, 2 * np.eye(4 * v))
    print(f"v={v}: complementary={ok} Hadamard={had} skew={skew} rows={[int(x.sum()) for x in seqs]}", flush=True)
    if ok and had and skew:
        with open(W.hits(), 'a') as f:
            f.write(json.dumps(dict(n=v, gens=[h], spin_multiplier=mu, pattern=[int(x.sum()) for x in seqs],
                                    hadamard=True, skew=True,
                                    seqs=[''.join('+' if x > 0 else '-' for x in s) for s in seqs])) + '\n')


if __name__ == '__main__':
    main()

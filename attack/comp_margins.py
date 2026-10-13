"""Margins (compressions to Z_p and Z_q) of a solution line "SOL X1 X2 ..." for comp_ls input files.

  python3 attack/comp_margins.py "SOL ..." p q types weights T > file      (types like s,f,f; weights like 1,2,1)
"""
import sys


def margins(seq, p, q):
    x = [1 if c == '+' else -1 for c in seq]
    v = len(x)
    rows = [sum(x[t] for t in range(v) if t % p == a) for a in range(p)]
    cols = [sum(x[t] for t in range(v) if t % q == b) for b in range(q)]
    return rows, cols


def main():
    line, p, q, types, weights, T = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4].split(','), sys.argv[5].split(','), int(sys.argv[6])
    seqs = line.split()[1:1 + len(types)]
    v = len(seqs[0])
    assert p * q == v
    out = [f"{v} {p} {q} {len(seqs)} {T}"]
    for s, t, w in zip(seqs, types, weights):
        r, c = margins(s, p, q)
        out.append(f"{t} {w} " + ' '.join(map(str, r)) + ' ' + ' '.join(map(str, c)))
    print('\n'.join(out))


if __name__ == '__main__':
    main()

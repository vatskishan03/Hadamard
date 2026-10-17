#!/usr/bin/env python3
"""T-sequences of length 5t from base sequences BS(m + 1, m), t = 2m + 1: Yang's multiplication by the Yang number 5
(Seberry and Yamada 1992, Theorem 5.23, with the Golay pair (+, +), (+, -) of length 2).

With E, F, G, H base sequences of lengths m + 1, m + 1, m, m put A = (E + F)/2, B = (E - F)/2, C = (G + H)/2,
D = (G - H)/2, 0 = zeros of length t split as (m + 1, m), and X* the reverse of X. Then

    X = (A, C; 0; A, -C; 0; B*, 0)      Y = (B, D; 0; B, -D; 0; -A*, 0)
    Z = (0; A, C; 0; -A, C; 0, D*)      W = (0; B, D; 0; -B, D; 0, -C*)

are T-sequences of length 5t (four ternary sequences with disjoint supports covering every position and zero
aperiodic autocorrelation sum off 0). The signs are those of the survey's template that work; the program checks the
result, and with --search it re-derives them by trying all 2^12 sign choices.

    python3 search/yang5.py 5 --turyn t-matrices/T55.txt   # m = 5: Turyn's BS(6, 5) -> T-matrices of order 55

Without --turyn the base sequences are those of search/paley_blocks.py; --turyn uses the Turyn sequences of length 6
(Turyn 1974, as stored in SageMath's turyn_sequences_smallcases), which give t-matrices/T55.txt.
"""
import itertools, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paley_blocks as PB


def is_ts(T):
    T = np.asarray(T, np.int64)
    n = T.shape[1]
    if not (np.abs(T).sum(0) == 1).all():
        return False
    return all(sum(int(r[:n - k] @ r[k:]) for r in T) == 0 for k in range(1, n))


TURYN = {5: ([1, 1, 1, -1, -1, -1], [1, 1, -1, 1, -1, 1], [1, 1, -1, 1, 1], [1, 1, -1, 1, 1])}


def yang5(m, signs=(1, 1, 1, -1, 1, 1, -1, 1, 1, -1, 1, -1), turyn=False):
    e, f, g, h = (np.asarray(x, np.int64) for x in (TURYN[m] if turyn else PB.base_sequences(m)))
    nac = lambda x, k: int(x[:len(x) - k] @ x[k:])
    assert all(sum(nac(x, k) for x in (e, f, g, h)) == 0 for k in range(1, m + 1)), "not base sequences"
    A, B, C, D = (e + f) // 2, (e - f) // 2, (g + h) // 2, (g - h) // 2
    Z1, Z0 = np.zeros(m + 1, np.int64), np.zeros(m, np.int64)
    e1, e2, e3, e4, d1, d2, d3, d4, s1, s2, s3, s4 = signs
    cat = np.concatenate
    X = cat([e1 * A, e2 * C, Z1, Z0, e3 * A, e4 * C, Z1, Z0, s1 * B[::-1], Z0])
    Y = cat([e1 * B, e2 * D, Z1, Z0, e3 * B, e4 * D, Z1, Z0, s2 * A[::-1], Z0])
    Zs = cat([Z1, Z0, d1 * A, d2 * C, Z1, Z0, d3 * A, d4 * C, Z1, s3 * D[::-1]])
    W = cat([Z1, Z0, d1 * B, d2 * D, Z1, Z0, d3 * B, d4 * D, Z1, s4 * C[::-1]])
    return np.array([X, Y, Zs, W])


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__); raise SystemExit(0)
    m = int(sys.argv[1])
    tur = "--turyn" in sys.argv
    args = [a for a in sys.argv[2:] if not a.startswith("--")]
    if "--search" in sys.argv:
        good = [s for s in itertools.product((1, -1), repeat=12) if is_ts(yang5(m, s, tur))]
        print(f"{len(good)} of 4096 sign choices give T-sequences of length {5 * (2 * m + 1)}")
    T = yang5(m, turyn=tur)
    assert is_ts(T), "not T-sequences"
    t = T.shape[1]
    paf = sum(np.array([int(r @ np.roll(r, -k)) for k in range(t)]) for r in T)
    assert paf[0] == t and not paf[1:].any()
    out = "\n".join("".join({1: "+", -1: "-", 0: "0"}[int(v)] for v in r) for r in T) + "\n"
    if args:
        Path(args[0]).write_text(out)
        print(f"wrote {args[0]}: T-matrices of order {t}")
    else:
        print(out, end="")

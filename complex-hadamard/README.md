# Complex Hadamard matrices of orders 2N, N odd

A quaternary (complex) Hadamard matrix of order n has entries 1, −1, i, −i and satisfies K K* = n I. Its order is 1
or even, and it is conjectured that one exists for every even order (Turyn). Orders divisible by 4 are covered by real
Hadamard matrices. The open cases are therefore the orders 2N with N odd.

Recent progress at the bottom of that range:
- Szöllősi constructed order 94 (arXiv:2603.09572, March 2026), noting that "the next outstanding order seems to be"
  118.
- Kharaghani, Mohammadian, Tayfeh-Rezaie and Zaitsev (arXiv:2605.08661v2, September 2026) gave the first examples of
  orders 118, 130, 134, 206 and 266, from near-Williamson matrices.

This directory adds orders that, to our knowledge, had no construction.

## Construction

Szöllősi, Theorem 4 of arXiv:2603.09572. Let A, B, C, D be circulant ±1 matrices of odd order N with

    A Aᵀ + B Bᵀ + C Cᵀ + D Dᵀ = 4N I,          A and B symmetric,

and let R be the back-diagonal identity. Then

    K = [ (A+B)/2 + i(A−B)/2           (C+DR)/2 + i(C−DR)/2 ]
        [ (Cᵀ+DR)/2 − i(Cᵀ−DR)/2      −(A+B)/2 + i(A−B)/2  ]

is a quaternary Hadamard matrix of order 2N.

**Proof sketch.** Put X = (A+B)/2 + i(A−B)/2 and Y = (C+DR)/2 + i(C−DR)/2. Their entries lie in {±1, ±i}. X is
symmetric and circulant, so it commutes with C, D and R, and K = [[X, Y], [Y*, −X*]].
- Since A, B and C, DR are amicable pairs, X X* + Y Y* = (A Aᵀ + B Bᵀ + C Cᵀ + D Dᵀ)/2 = 2N I.
- X and Y are normal, so X* X + Y* Y = 2N I as well.
- XY = YX makes the two block rows orthogonal.

Williamson and near-Williamson matrices are the special cases with three or four symmetric blocks. Allowing two
arbitrary circulant blocks enlarges the search space by far.

## Search

All four blocks are constant on the orbits of a multiplier group H ≤ ℤ_N* with −1 ∉ H, so C and D are H-invariant.
A and B are symmetric, hence constant on the orbits of ⟨H, −1⟩. The condition Σ PAF = 0 off zero then holds on the
⟨H, −1⟩-classes of ℤ_N \ {0}.

For each row-sum pattern (a, b, c, d) with a² + b² + c² + d² = 4N:
- the four blocks are listed completely (`attack/bnw_enum.c`);
- the lists are joined exactly (`attack/join4x.c`, pairing each symmetric list with an arbitrary one);
- every hit is turned into K and checked exactly, and the real Goethals–Seidel matrix of order 4N is checked as well.

Code: `attack/ssxx.py` (triage, search), `attack/ssxx_sweep.py` (the candidate orders and the sweep),
`attack/ssxx_certify.py` (certificates).

## Matrices

| order 2N | N | H | row sums (a, b; c, d) | solutions in the structure (exact join) | file |
|---:|---:|---|---|---:|---|
| 142 | 71 | ⟨5⟩, order 5 | 9, 9; 1, 11 | 5992 (and 6076 with 9, 11; 1, 9) | [C142.json](C142.json) |
| 154 | 77 = 7 · 11 | ⟨15⟩, order 5 (15 ≡ 1 mod 7, the squares mod 11) | 3, 3; 1, 17 | 2352 (also 2352, 2352 and 1176 for three other patterns; the Gaussian model predicts 10^−1.7) | [C154.json](C154.json) |
| 222 | 111 = 3 · 37 | from the near-Williamson matrices of order 111 in [`near-williamson/NW111.txt`](../near-williamson/NW111.txt) (order-3 twist, `attack/wt_twist.py`), with (B, C) as the symmetric pair and (A, D) as the other | 1, 1; 21, 1 (in the order B, C; A, D) | 36 near-Williamson quadruples with these rows | [C222.json](C222.json) |

Each certificate gives N, the multiplier group, the row sums, the first rows of A, B, C, D (`+`, `-`) and the
SHA-256 of the canonical text of K: 2N lines of 2N digits, digit e standing for the entry iᵉ.

```
python3 verify/verify_complex.py            # all certificates
python3 verify/verify_complex.py --write 142 # also writes complex-hadamard/C142.txt
```

This rebuilds K from the four sequences alone and checks every entry and K K* = 2N I exactly, in integer arithmetic
on the real and imaginary parts.

## Status

To our knowledge, these orders had no construction. Our reading of the literature counts 2N as known when one of the
following gives it:

- **Williamson-type matrices of order N** (Turyn; construction (3) of arXiv:2605.08661). Sources:
  - Table A.1 of Seberry and Yamada (1992), parsed in [`attack/data/wt_known_1992.txt`](../attack/data/wt_known_1992.txt);
  - Đoković's good matrices of orders 35 and 127 (1993);
  - the 2026 near-Williamson orders: every odd N ≤ 69, and 73, 79, 91, 93, 95, 97, 103, 109, 127, 133, 181.
- **A Paley conference matrix of order 2N** (2N − 1 a prime power), via C + iI.
- **Complex Golay sequences of length N** (N = 1, 3, 5, 11, 13).
- **A quaternary Legendre pair of length N − 1:**
  - every even length up to 100 except 70, 76, 88, 92, 94 and 100 (Kotsireas–Winterhof; Kotsireas–Koutschan–Winterhof;
    Lebedev and Verma 2026);
  - the two families of Jedwab and Pender (2025).
- **Szöllősi's order 94.**

Orders 142 and 154 are the two smallest orders that none of these give.
- 141 = 3 · 47 and 153 = 3² · 17 are not prime powers.
- 71 and 77 are absent from Table A.1 as orders of four Williamson-type matrices (71 is blank; 77 has only
  Williamson-type matrices of order 154), and from the 2026 near-Williamson list.
- Quaternary Legendre pairs of lengths 70 and 76 are open; Lebedev's exhaustive searches of the reversible family
  at these lengths found none.

The 2026 papers above cite Table 11.2 of Seberry and Yamada's 2020 book for the list of open orders. We could not
consult it, so these claims carry that caveat.

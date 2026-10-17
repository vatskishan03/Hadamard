# Hadamard matrix of order 26740

26740 = 4 · 6685, with 6685 = 5 · 7 · 191.

## Construction

Paley's skew-Hadamard matrix of order 5348 (q = 5347, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 1337 once its points are listed along the four orbits of a torus element of PSL(2, 5347) and rescaled by
signs. Multiplying it by T-matrices of order 5 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 6685 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 5 (from base sequences BS(3, 2), built by the program): first rows `++000`, `00+00`, `000+-`,
  `00000`.
- Row sums of the four sequences: 5, −5, −67, 149 (5² + (−5)² + (−67)² + 149² = 26740).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 5347 5` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 6685 has t = 3: the smallest known Hadamard matrix with
odd part 6685 had order 2³ · 6685. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 5348 is equivalent
  to an array of sixteen circulant blocks of order 1337 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and
  Yamada (1990, Theorem 4) multiply such an array by T-matrices of order 5.
- Paley's constructions do not apply: 4n − 1 = 26739 = 3² · 2971 and 2n − 1 = 13369 = 29 · 461 are not prime powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 6686 cannot exist (6685 is not a sum
  of two squares); Yamada–Kiyasu does not apply: n − 2 = 6683 = 41 · 163 is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(26740, existence=True)`
  returns Unknown.

We found no construction or statement of order 26740. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 6685 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 26740
```

SHA-256 of `H26740.csv`: `465a61bf711406be479e1240207bf530803a6ce861b6bfb61ea8503722ca398d`

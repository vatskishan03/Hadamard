# Hadamard matrix of order 15764

15764 = 4 · 3941, with 3941 = 7 · 563.

## Construction

Paley's skew-Hadamard matrix of order 2252 (q = 2251, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 563 once its points are listed along the four orbits of a torus element of PSL(2, 2251) and rescaled by
signs. Multiplying it by T-matrices of order 7 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 3941 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = −1.
- T-matrices of order 7 (from base sequences BS(4, 3), built by the program): first rows `+000000`, `0++-000`,
  `0000+0+`, `00000+0`.
- Row sums of the four sequences: 27, 59, −77, 75 (27² + 59² + (−77)² + 75² = 15764).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 2251 7` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 3941 has t = 3: the smallest known Hadamard matrix with
odd part 3941 had order 2³ · 3941. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 2252 is equivalent
  to an array of sixteen circulant blocks of order 563 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 7.
- Paley's constructions do not apply: 4n − 1 = 15763 = 11 · 1433 and 2n − 1 = 7881 = 3 · 37 · 71 are not prime
  powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 3942 cannot exist (3941 is not a sum
  of two squares); Yamada–Kiyasu does not apply: n − 2 = 3939 = 3 · 13 · 101 is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(15764, existence=True)`
  returns Unknown.

We found no construction or statement of order 15764. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 3941 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 15764
```

SHA-256 of `H15764.csv`: `fc077676f5ab7c34374ee2bffc8d43f933049d05b80fae724ec40e539404895c`

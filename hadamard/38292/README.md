# Hadamard matrix of order 38292

38292 = 4 · 9573, with 9573 = 3 · 3191.

## Construction

Paley's skew-Hadamard matrix of order 12764 (q = 12763, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 3191 once its points are listed along the four orbits of a torus element of PSL(2, 12763) and rescaled by
signs. Multiplying it by T-matrices of order 3 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 9573 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 3 (from base sequences BS(2, 1), built by the program): first rows `+00`, `0+0`, `00+`, `000`.
- Row sums of the four sequences: 89, 43, −129, 109 (89² + 43² + (−129)² + 109² = 38292).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 12763 3` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 9573 has t = 3: the smallest known Hadamard matrix with
odd part 9573 had order 2³ · 9573. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 12764 is equivalent
  to an array of sixteen circulant blocks of order 3191 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and
  Yamada (1990, Theorem 4) multiply such an array by T-matrices of order 3.
- Paley's constructions do not apply: 4n − 1 = 38291 = 11 · 59² and 2n − 1 = 19145 = 5 · 7 · 547 are not prime
  powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 9574 cannot exist (9573 is not a sum
  of two squares); Yamada–Kiyasu does not apply: n − 2 = 9571 = 17 · 563 is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(38292, existence=True)`
  returns Unknown.

We found no construction or statement of order 38292. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 9573 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 38292
```

SHA-256 of `H38292.csv`: `87c8e9ad773281aeaeeca473e5f68706961b4e42c928ad4799c04bac4f37f9bf`

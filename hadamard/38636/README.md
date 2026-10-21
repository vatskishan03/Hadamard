# Hadamard matrix of order 38636

38636 = 4 · 9659, with 9659 = 13 · 743.

## Construction

Paley's skew-Hadamard matrix of order 2972 (q = 2971, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 743 once its points are listed along the four orbits of a torus element of PSL(2, 2971) and rescaled by
signs. Multiplying it by T-matrices of order 13 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 9659 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = +1.
- T-matrices of order 13 (from base sequences BS(7, 6), built by the program): first rows `+0+0+00000000`,
  `0+0+0-+000000`, `0000000+0000-`, `00000000++--0`.
- Row sums of the four sequences: 21, −25, 133, −141 (21² + (−25)² + 133² + (−141)² = 38636).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 2971 13` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 9659 has t = 3: the smallest known Hadamard matrix with
odd part 9659 had order 2³ · 9659. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 2972 is equivalent
  to an array of sixteen circulant blocks of order 743 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 13.
- Paley's constructions do not apply: 4n − 1 = 38635 = 5 · 7727 and 2n − 1 = 19317 = 3 · 47 · 137 are not prime
  powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4); Yamada–Kiyasu does not apply: n − 2 = 9657 = 3² · 29 · 37 is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(38636, existence=True)`
  returns Unknown.

We found no construction or statement of order 38636. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 9659 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 38636
```

SHA-256 of `H38636.csv`: `8f64774a3a1cfd3e6e2bf4fba4a7322a32b5d6089b81033db772f95009f9f82e`

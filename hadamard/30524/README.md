# Hadamard matrix of order 30524

30524 = 4 · 7631, with 7631 = 13 · 587.

## Construction

Paley's skew-Hadamard matrix of order 2348 (q = 2347, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 587 once its points are listed along the four orbits of a torus element of PSL(2, 2347) and rescaled by
signs. Multiplying it by T-matrices of order 13 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 7631 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 13 (from base sequences BS(7, 6), built by the program): first rows `+0+0+00000000`,
  `0+0+0-+000000`, `0000000+0000-`, `00000000++--0`.
- Row sums of the four sequences: −63, 101, 125, 27 ((−63)² + 101² + 125² + 27² = 30524).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 2347 13` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 7631 has t = 3: the smallest known Hadamard matrix with
odd part 7631 had order 2³ · 7631. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 2348 is equivalent
  to an array of sixteen circulant blocks of order 587 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 13.
- Paley's constructions do not apply: 4n − 1 = 30523 = 131 · 233 and 2n − 1 = 15261 = 3 · 5087 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4); Yamada–Kiyasu does not apply: n − 2 = 7629 = 3 · 2543 is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(30524, existence=True)`
  returns Unknown.

We found no construction or statement of order 30524. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 7631 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 30524
```

SHA-256 of `H30524.csv`: `82bfe236ad3f062430e1db943b75d2cc55e8b9125bfed0a9410092ca2d302aed`

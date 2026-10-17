# Hadamard matrix of order 24836

24836 = 4 · 6209, with 6209 = 7 · 887.

## Construction

Paley's skew-Hadamard matrix of order 3548 (q = 3547, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 887 once its points are listed along the four orbits of a torus element of PSL(2, 3547) and rescaled by
signs. Multiplying it by T-matrices of order 7 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 6209 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 7 (from base sequences BS(4, 3), built by the program): first rows `+000000`, `0++-000`,
  `0000+0+`, `00000+0`.
- Row sums of the four sequences: 71, 75, −101, 63 (71² + 75² + (−101)² + 63² = 24836).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 3547 7` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 6209 has t = 8: the smallest known Hadamard matrix with
odd part 6209 had order 2⁸ · 6209. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 3548 is equivalent
  to an array of sixteen circulant blocks of order 887 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 7.
- Paley's constructions do not apply: 4n − 1 = 24835 = 5 · 4967 and 2n − 1 = 12417 = 3 · 4139 are not prime powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 6210 cannot exist (6209 is not a sum
  of two squares); Yamada–Kiyasu does not apply: n − 2 = 6207 = 3 · 2069 is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(24836, existence=True)`
  returns Unknown.

We found no construction or statement of order 24836. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 6209 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 24836
```

SHA-256 of `H24836.csv`: `f660c3b9dc0aaffaae7abc1b90cf20d4d78b905bebe7cb72143f375e950ceaf7`

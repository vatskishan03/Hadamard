# Hadamard matrix of order 25860

25860 = 4 · 6465, with 6465 = 3 · 5 · 431.

## Construction

Paley's skew-Hadamard matrix of order 5172 (q = 5171, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 1293 once its points are listed along the four orbits of a torus element of PSL(2, 5171) and rescaled by
signs. Multiplying it by T-matrices of order 5 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 6465 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 6; sign normalisation c = −1.
- T-matrices of order 5 (from base sequences BS(3, 2), built by the program): first rows `++000`, `00+00`, `000+-`,
  `00000`.
- Row sums of the four sequences: 33, −61, 83, −119 (33² + (−61)² + 83² + (−119)² = 25860).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 5171 5` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 6465 has t = 3: the smallest known Hadamard matrix with
odd part 6465 had order 2³ · 6465. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 5172 is equivalent
  to an array of sixteen circulant blocks of order 1293 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and
  Yamada (1990, Theorem 4) multiply such an array by T-matrices of order 5.
- Paley's constructions do not apply: 4n − 1 = 25859 = 19 · 1361 and 2n − 1 = 12929 = 7 · 1847 are not prime powers.
- Miyamoto and Yamada–Kiyasu do not apply: 6465 is not a prime power, and neither is n − 2 = 6463 = 23 · 281.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(25860, existence=True)`
  returns Unknown.

We found no construction or statement of order 25860. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 6465 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 25860
```

SHA-256 of `H25860.csv`: `c58af96eb30b55ddd4db04ead3e4c42c5e5228f125599e10c063c6bb1fecffcd`

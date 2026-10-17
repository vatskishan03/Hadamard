# Hadamard matrix of order 26124

26124 = 4 · 6531, with 6531 = 3 · 7 · 311.

## Construction

Paley's skew-Hadamard matrix of order 8708 (q = 8707, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 2177 once its points are listed along the four orbits of a torus element of PSL(2, 8707) and rescaled by
signs. Multiplying it by T-matrices of order 3 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 6531 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 3 (from base sequences BS(2, 1), built by the program): first rows `+00`, `0+0`, `00+`, `000`.
- Row sums of the four sequences: 105, −113, 11, −47 (105² + (−113)² + 11² + (−47)² = 26124).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 8707 3` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 6531 has t = 5: the smallest known Hadamard matrix with
odd part 6531 had order 2⁵ · 6531. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 8708 is equivalent
  to an array of sixteen circulant blocks of order 2177 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and
  Yamada (1990, Theorem 4) multiply such an array by T-matrices of order 3.
- Paley's constructions do not apply: 4n − 1 = 26123 = 151 · 173 and 2n − 1 = 13061 = 37 · 353 are not prime powers.
- Miyamoto does not apply (6531 is not a prime power). Yamada's theorem with q = n − 2 = 6529 needs a symmetric
  conference matrix of order 3266; 3265 = 5 · 653 is not a prime power and we know no such matrix.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(26124, existence=True)`
  returns Unknown.

We found no construction or statement of order 26124. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 6531 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 26124
```

SHA-256 of `H26124.csv`: `0f744ab7c081098a79bf44a0771e381ea91adc2ed293423d349113bb8b8788ec`

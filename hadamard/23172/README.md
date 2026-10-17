# Hadamard matrix of order 23172

23172 = 4 · 5793, with 5793 = 3 · 1931.

## Construction

Paley's skew-Hadamard matrix of order 7724 (q = 7723, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 1931 once its points are listed along the four orbits of a torus element of PSL(2, 7723) and rescaled by
signs. Multiplying it by T-matrices of order 3 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 5793 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 3 (from base sequences BS(2, 1), built by the program): first rows `+00`, `0+0`, `00+`, `000`.
- Row sums of the four sequences: 1, 93, −91, 79 (1² + 93² + (−91)² + 79² = 23172).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 7723 3` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 5793 has t = 3: the smallest known Hadamard matrix with
odd part 5793 had order 2³ · 5793. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 7724 is equivalent
  to an array of sixteen circulant blocks of order 1931 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and
  Yamada (1990, Theorem 4) multiply such an array by T-matrices of order 3.
- Paley's constructions do not apply: 4n − 1 = 23171 = 17 · 29 · 47 and 2n − 1 = 11585 = 5 · 7 · 331 are not prime
  powers.
- Miyamoto and Yamada–Kiyasu do not apply: 5793 is not a prime power, and n − 2 = 5791 is ≡ 7 (mod 8).
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(23172, existence=True)`
  returns Unknown.

We found no construction or statement of order 23172. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 5793 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 23172
```

SHA-256 of `H23172.csv`: `da08784e50649e1486f384f1ac00b54b78435424e157f12b26aad453736c3032`

# Hadamard matrix of order 35484

35484 = 4 · 8871, with 8871 = 3 · 2957.

## Construction

Paley's skew-Hadamard matrix of order 11828 (q = 11827, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 2957 once its points are listed along the four orbits of a torus element of PSL(2, 11827) and rescaled by
signs. Multiplying it by T-matrices of order 3 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 8871 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 3 (from base sequences BS(2, 1), built by the program): first rows `+00`, `0+0`, `00+`, `000`.
- Row sums of the four sequences: 1, 131, −129, 41 (1² + 131² + (−129)² + 41² = 35484).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 11827 3` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 8871 has t = 3: the smallest known Hadamard matrix with
odd part 8871 had order 2³ · 8871. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 11828 is equivalent
  to an array of sixteen circulant blocks of order 2957 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and
  Yamada (1990, Theorem 4) multiply such an array by T-matrices of order 3.
- Paley's constructions do not apply: 4n − 1 = 35483 = 7 · 37 · 137 and 2n − 1 = 17741 = 113 · 157 are not prime
  powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4); Yamada–Kiyasu does not apply: n − 2 = 8869 = 7² · 181 is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(35484, existence=True)`
  returns Unknown.

We found no construction or statement of order 35484. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 8871 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 35484
```

SHA-256 of `H35484.csv`: `d5b72faf3b23a6642fa5286a0874f8d85146288211bf7a34496408667676d0df`

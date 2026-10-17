# Hadamard matrix of order 38228

38228 = 4 · 9557, with 9557 = 19 · 503.

## Construction

Paley's skew-Hadamard matrix of order 2012 (q = 2011, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 503 once its points are listed along the four orbits of a torus element of PSL(2, 2011) and rescaled by
signs. Multiplying it by T-matrices of order 19 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 9557 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = +1.
- T-matrices of order 19 ([`t-matrices/T19.txt`](../../t-matrices/T19.txt), found by `search/tmat_mult.py 19 7`):
  first rows `++00000+-00+-00000-`, `0000000000000000000`, `00++0000000000+0000`, `0000+++00+-00-0-++0`.
- Row sums of the four sequences: 1, −97, −117, 123 (1² + (−97)² + (−117)² + 123² = 38228).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 2011 19
t-matrices/T19.txt` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 9557 has t = 3: the smallest known Hadamard matrix with
odd part 9557 had order 2³ · 9557. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 2012 is equivalent
  to an array of sixteen circulant blocks of order 503 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 19.
- Paley's constructions do not apply: 4n − 1 = 38227 = 7 · 43 · 127 and 2n − 1 = 19113 = 3 · 23 · 277 are not prime
  powers.
- Miyamoto and Yamada–Kiyasu do not apply: 9557 is not a prime power, and neither is n − 2 = 9555 = 3 · 5 · 7² · 13.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(38228, existence=True)`
  returns Unknown.

We found no construction or statement of order 38228. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 9557 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 38228
```

SHA-256 of `H38228.csv`: `82a111660ed2db5edd2760ec859ff07caa5eec4ada1a5eb1c1b60cd49b4db601`

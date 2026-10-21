# Hadamard matrix of order 28932

28932 = 4 · 7233, with 7233 = 3 · 2411.

## Construction

Paley's skew-Hadamard matrix of order 9644 (q = 9643, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 2411 once its points are listed along the four orbits of a torus element of PSL(2, 9643) and rescaled by
signs. Multiplying it by T-matrices of order 3 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 7233 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 3 (from base sequences BS(2, 1), built by the program): first rows `+00`, `0+0`, `00+`, `000`.
- Row sums of the four sequences: −57, 1, 59, 149 ((−57)² + 1² + 59² + 149² = 28932).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 9643 3` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 7233 has t = 3: the smallest known Hadamard matrix with
odd part 7233 had order 2³ · 7233. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 9644 is equivalent
  to an array of sixteen circulant blocks of order 2411 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and
  Yamada (1990, Theorem 4) multiply such an array by T-matrices of order 3.
- Paley's constructions do not apply: 4n − 1 = 28931 = 7 · 4133 and 2n − 1 = 14465 = 5 · 11 · 263 are not prime
  powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 7234 cannot exist (7233 is not a sum
  of two squares); Yamada–Kiyasu does not apply: n − 2 = 7231 = 7 · 1033 is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(28932, existence=True)`
  returns Unknown.

We found no construction or statement of order 28932. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 7233 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 28932
```

SHA-256 of `H28932.csv`: `44e41a91c92dcc8a2addb05b33f7ed487ae535dc45303e8e04df2109b8e50665`

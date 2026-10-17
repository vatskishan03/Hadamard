# Hadamard matrix of order 19284

19284 = 4 · 4821, with 4821 = 3 · 1607.

## Construction

Paley's skew-Hadamard matrix of order 6428 (q = 6427, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 1607 once its points are listed along the four orbits of a torus element of PSL(2, 6427) and rescaled by
signs. Multiplying it by T-matrices of order 3 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 4821 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 3 (from base sequences BS(2, 1), built by the program): first rows `+00`, `0+0`, `00+`, `000`.
- Row sums of the four sequences: −15, 103, −85, 35 ((−15)² + 103² + (−85)² + 35² = 19284).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 6427 3` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 4821 has t = 3: the smallest known Hadamard matrix with
odd part 4821 had order 2³ · 4821. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 6428 is equivalent
  to an array of sixteen circulant blocks of order 1607 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and
  Yamada (1990, Theorem 4) multiply such an array by T-matrices of order 3.
- Paley's constructions do not apply: 4n − 1 = 19283 = 11 · 1753 and 2n − 1 = 9641 = 31 · 311 are not prime powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 4822 cannot exist (4821 is not a sum
  of two squares); Yamada–Kiyasu does not apply: n − 2 = 4819 = 61 · 79 is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(19284, existence=True)`
  returns Unknown.

We found no construction or statement of order 19284. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 4821 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 19284
```

SHA-256 of `H19284.csv`: `bc2d3e281cb3a5706f7fcf7ea28b13aac343f1a4af0baca242dba75096aca787`

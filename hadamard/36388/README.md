# Hadamard matrix of order 36388

36388 = 4 · 9097, with 9097 = 11 · 827.

## Construction

Paley's skew-Hadamard matrix of order 3308 (q = 3307, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 827 once its points are listed along the four orbits of a torus element of PSL(2, 3307) and rescaled by
signs. Multiplying it by T-matrices of order 11 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 9097 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 11 ([`t-matrices/T11.txt`](../../t-matrices/T11.txt), T-sequences of length 11 from a Golay
  pair of length 10 (Turyn 1974)): first rows `+0000000000`, `0++00+0+-00`, `000++0-00-+`, `00000000000`.
- Row sums of the four sequences: 109, 19, −155, 11 (109² + 19² + (−155)² + 11² = 36388).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 3307 11
t-matrices/T11.txt` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 9097 has t = 3: the smallest known Hadamard matrix with
odd part 9097 had order 2³ · 9097. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 3308 is equivalent
  to an array of sixteen circulant blocks of order 827 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 11.
- Paley's constructions do not apply: 4n − 1 = 36387 = 3² · 13 · 311 and 2n − 1 = 18193 = 7 · 23 · 113 are not prime
  powers.
- Miyamoto and Yamada–Kiyasu do not apply: 9097 is not a prime power, and neither is n − 2 = 9095 = 5 · 17 · 107.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(36388, existence=True)`
  returns Unknown.

We found no construction or statement of order 36388. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 9097 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 36388
```

SHA-256 of `H36388.csv`: `eb3518dd7fbe40b05aa30fe79348d94a07ce725fcc3665e108a7ee71a11437b0`

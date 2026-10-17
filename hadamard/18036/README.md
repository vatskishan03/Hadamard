# Hadamard matrix of order 18036

18036 = 4 · 4509, with 4509 = 3³ · 167.

## Construction

Paley's skew-Hadamard matrix of order 6012 (q = 6011, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 1503 once its points are listed along the four orbits of a torus element of PSL(2, 6011) and rescaled by
signs. Multiplying it by T-matrices of order 3 (Seberry–Yamada 1990, Theorem 4) gives four ±1 sequences on the group
ℤ₃ × ℤ₁₅₀₃ of order 4509 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The
construction and its proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 6; sign normalisation c = −1.
- T-matrices of order 3 (from base sequences BS(2, 1), built by the program): first rows `+00`, `0+0`, `00+`, `000`.
- Row sums of the four sequences: −49, −43, 95, −69 ((−49)² + (−43)² + 95² + (−69)² = 18036).
- Since gcd(s, t) = 3, the sequences live on ℤ_3 × ℤ_1503, not on ℤ_4509. `sequences.txt` lists the entry at (i, a)
  in position i · 1503 + a, and the Goethals–Seidel array uses the group-developed matrices and the reversal x ↦ −x
  of the group.

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 6011 3` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 4509 has t = 3: the smallest known Hadamard matrix with
odd part 4509 had order 2³ · 4509. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 6012 is equivalent
  to an array of sixteen circulant blocks of order 1503 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and
  Yamada (1990, Theorem 4) multiply such an array by T-matrices of order 3.
- Paley's constructions do not apply: 4n − 1 = 18035 = 5 · 3607 and 2n − 1 = 9017 = 71 · 127 are not prime powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 4510 cannot exist (4509 is not a sum
  of two squares); Yamada–Kiyasu does not apply: n − 2 = 4507 is ≡ 3 (mod 8).
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(18036, existence=True)`
  returns Unknown.

We found no construction or statement of order 18036. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 4509 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 18036
```

SHA-256 of `H18036.csv`: `eef6724b13e56e95c5ab00d58066c99266c3d544a568b5e9f9db3732fa4a6d8e`

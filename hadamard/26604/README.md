# Hadamard matrix of order 26604

26604 = 4 · 6651, with 6651 = 3² · 739.

## Construction

Paley's skew-Hadamard matrix of order 8868 (q = 8867, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 2217 once its points are listed along the four orbits of a torus element of PSL(2, 8867) and rescaled by
signs. Multiplying it by T-matrices of order 3 (Seberry–Yamada 1990, Theorem 4) gives four ±1 sequences on the group
ℤ₃ × ℤ₂₂₁₇ of order 6651 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The
construction and its proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 3 (from base sequences BS(2, 1), built by the program): first rows `+00`, `0+0`, `00+`, `000`.
- Row sums of the four sequences: 111, 1, −109, 49 (111² + 1² + (−109)² + 49² = 26604).
- Since gcd(s, t) = 3, the sequences live on ℤ_3 × ℤ_2217, not on ℤ_6651. `sequences.txt` lists the entry at (i, a)
  in position i · 2217 + a, and the Goethals–Seidel array uses the group-developed matrices and the reversal x ↦ −x
  of the group.

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 8867 3` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 6651 has t = 7: the smallest known Hadamard matrix with
odd part 6651 had order 2⁷ · 6651. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 8868 is equivalent
  to an array of sixteen circulant blocks of order 2217 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and
  Yamada (1990, Theorem 4) multiply such an array by T-matrices of order 3.
- Paley's constructions do not apply: 4n − 1 = 26603 = 37 · 719 and 2n − 1 = 13301 = 47 · 283 are not prime powers.
- Miyamoto and Yamada–Kiyasu do not apply: 6651 is not a prime power, and neither is n − 2 = 6649 = 61 · 109.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(26604, existence=True)`
  returns Unknown.

We found no construction or statement of order 26604. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 6651 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 26604
```

SHA-256 of `H26604.csv`: `8b9101c1b15bbd088664d5dfd98eef424dba06cdc6c1a6d0c41df045934e2c94`

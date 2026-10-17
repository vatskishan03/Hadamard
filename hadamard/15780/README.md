# Hadamard matrix of order 15780

15780 = 4 · 3945, with 3945 = 3 · 5 · 263.

## Construction

Paley's skew-Hadamard matrix of order 1052 (q = 1051, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 263 once its points are listed along the four orbits of a torus element of PSL(2, 1051) and rescaled by
signs. Multiplying it by T-matrices of order 15 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 3945 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = −1.
- T-matrices of order 15 ([`t-matrices/T15.txt`](../../t-matrices/T15.txt), T-sequences of length 15 as stored in
  SageMath's `t_sequences` module): first rows `++0+-0-00000000`, `00-00+0-0000000`, `00000000+00000+`,
  `000000000++-++0`.
- Row sums of the four sequences: −117, −35, 5, 29 ((−117)² + (−35)² + 5² + 29² = 15780).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 1051 15
t-matrices/T15.txt` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 3945 has t = 3: the smallest known Hadamard matrix with
odd part 3945 had order 2³ · 3945. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 1052 is equivalent
  to an array of sixteen circulant blocks of order 263 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 15.
- Paley's constructions do not apply: 4n − 1 = 15779 = 31 · 509 and 2n − 1 = 7889 = 7³ · 23 are not prime powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 3946 cannot exist (3945 is not a sum
  of two squares); Yamada–Kiyasu does not apply: n − 2 = 3943 is ≡ 7 (mod 8).
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(15780, existence=True)`
  returns Unknown.

We found no construction or statement of order 15780. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 3945 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 15780
```

SHA-256 of `H15780.csv`: `73b3f36078c85898b1546101ea4ebb8bedb0b07cca50b57c2316db38312dc312`

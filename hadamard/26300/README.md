# Hadamard matrix of order 26300

26300 = 4 · 6575, with 6575 = 5² · 263.

## Construction

Paley's skew-Hadamard matrix of order 1052 (q = 1051, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 263 once its points are listed along the four orbits of a torus element of PSL(2, 1051) and rescaled by
signs. Multiplying it by T-matrices of order 25 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 6575 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = −1.
- T-matrices of order 25 ([`t-matrices/T25.txt`](../../t-matrices/T25.txt), T-sequences of length 25 as stored in
  SageMath's `t_sequences` module): first rows `+++0-+-+-0++0000000000000`, `000+00000+00+000000000000`,
  `0000000000000+++-0+-0+---`, `00000000000000000+00-0000`.
- Row sums of the four sequences: 43, −49, 21, 147 (43² + (−49)² + 21² + 147² = 26300).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 1051 25
t-matrices/T25.txt` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 6575 has t = 4: the smallest known Hadamard matrix with
odd part 6575 had order 2⁴ · 6575. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 1052 is equivalent
  to an array of sixteen circulant blocks of order 263 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 25.
- Paley's constructions do not apply: 4n − 1 = 26299 = 7 · 13 · 17² and 2n − 1 = 13149 = 3³ · 487 are not prime
  powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4); Yamada–Kiyasu does not apply: n − 2 = 6573 = 3 · 7 · 313 is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(26300, existence=True)`
  returns Unknown.

We found no construction or statement of order 26300. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 6575 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 26300
```

SHA-256 of `H26300.csv`: `0c0513e2c71daac6b31af1b15c3bf058df1c86ebfe6c0bfa540978f7877be45f`

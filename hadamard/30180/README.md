# Hadamard matrix of order 30180

30180 = 4 · 7545, with 7545 = 3 · 5 · 503.

## Construction

Paley's skew-Hadamard matrix of order 2012 (q = 2011, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 503 once its points are listed along the four orbits of a torus element of PSL(2, 2011) and rescaled by
signs. Multiplying it by T-matrices of order 15 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 7545 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = +1.
- T-matrices of order 15 ([`t-matrices/T15.txt`](../../t-matrices/T15.txt), T-sequences of length 15 as stored in
  SageMath's `t_sequences` module): first rows `++0+-0-00000000`, `00-00+0-0000000`, `00000000+00000+`,
  `000000000++-++0`.
- Row sums of the four sequences: −33, −89, −127, 71 ((−33)² + (−89)² + (−127)² + 71² = 30180).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 2011 15
t-matrices/T15.txt` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 7545 has t = 3: the smallest known Hadamard matrix with
odd part 7545 had order 2³ · 7545. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 2012 is equivalent
  to an array of sixteen circulant blocks of order 503 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 15.
- Paley's constructions do not apply: 4n − 1 = 30179 = 103 · 293 and 2n − 1 = 15089 = 79 · 191 are not prime powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 7546 cannot exist (7545 is not a sum
  of two squares); Yamada–Kiyasu does not apply: n − 2 = 7543 = 19 · 397 is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(30180, existence=True)`
  returns Unknown.

We found no construction or statement of order 30180. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 7545 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 30180
```

SHA-256 of `H30180.csv`: `b81add63c7520ff9a7bba3e4dc5b737b329d1aa2bf2b392bbf5adb5b62cdbacc`

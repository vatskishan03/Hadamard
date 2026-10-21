# Hadamard matrix of order 29108

29108 = 4 · 7277, with 7277 = 19 · 383.

## Construction

Paley's skew-Hadamard matrix of order 1532 (q = 1531, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 383 once its points are listed along the four orbits of a torus element of PSL(2, 1531) and rescaled by
signs. Multiplying it by T-matrices of order 19 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 7277 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = −1.
- T-matrices of order 19 ([`t-matrices/T19.txt`](../../t-matrices/T19.txt), found by `search/tmat_mult.py 19 7`):
  first rows `++00000+-00+-00000-`, `0000000000000000000`, `00++0000000000+0000`, `0000+++00+-00-0-++0`.
- Row sums of the four sequences: 55, −37, −117, 105 (55² + (−37)² + (−117)² + 105² = 29108).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 1531 19
t-matrices/T19.txt` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 7277 has t = 3: the smallest known Hadamard matrix with
odd part 7277 had order 2³ · 7277. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 1532 is equivalent
  to an array of sixteen circulant blocks of order 383 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 19.
- Paley's constructions do not apply: 4n − 1 = 29107 = 13 · 2239 and 2n − 1 = 14553 = 3³ · 7² · 11 are not prime
  powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 7278 cannot exist (7277 is not a sum
  of two squares); Yamada–Kiyasu does not apply: n − 2 = 7275 = 3 · 5² · 97 is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(29108, existence=True)`
  returns Unknown.

We found no construction or statement of order 29108. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 7277 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 29108
```

SHA-256 of `H29108.csv`: `39cbef22bfe4ae3785a840cd0960703052b00f9d5cab8c885883b971f2913ee4`

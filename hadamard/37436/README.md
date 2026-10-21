# Hadamard matrix of order 37436

37436 = 4 · 9359, with 9359 = 7² · 191.

## Construction

Paley's skew-Hadamard matrix of order 5348 (q = 5347, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 1337 once its points are listed along the four orbits of a torus element of PSL(2, 5347) and rescaled by
signs. Multiplying it by T-matrices of order 7 (Seberry–Yamada 1990, Theorem 4) gives four ±1 sequences on the group
ℤ₇ × ℤ₁₃₃₇ of order 9359 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The
construction and its proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 7 (from base sequences BS(4, 3), built by the program): first rows `+000000`, `0++-000`,
  `0000+0+`, `00000+0`.
- Row sums of the four sequences: −75, 141, −71, 83 ((−75)² + 141² + (−71)² + 83² = 37436).
- Since gcd(s, t) = 7, the sequences live on ℤ_7 × ℤ_1337, not on ℤ_9359. `sequences.txt` lists the entry at (i, a)
  in position i · 1337 + a, and the Goethals–Seidel array uses the group-developed matrices and the reversal x ↦ −x
  of the group.

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 5347 7` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 9359 has t = 3: the smallest known Hadamard matrix with
odd part 9359 had order 2³ · 9359. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 5348 is equivalent
  to an array of sixteen circulant blocks of order 1337 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and
  Yamada (1990, Theorem 4) multiply such an array by T-matrices of order 7.
- Paley's constructions do not apply: 4n − 1 = 37435 = 5 · 7487 and 2n − 1 = 18717 = 3 · 17 · 367 are not prime
  powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4); Yamada–Kiyasu does not apply: n − 2 = 9357 = 3 · 3119 is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(37436, existence=True)`
  returns Unknown.

We found no construction or statement of order 37436. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 9359 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 37436
```

SHA-256 of `H37436.csv`: `e67b77167c828c726825b14c5c754b6dc5430c1097ff8626a64bc5f61e23df6c`

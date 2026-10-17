# Hadamard matrix of order 17532

17532 = 4 · 4383, with 4383 = 3² · 487.

## Construction

Paley's skew-Hadamard matrix of order 5844 (q = 5843, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 1461 once its points are listed along the four orbits of a torus element of PSL(2, 5843) and rescaled by
signs. Multiplying it by T-matrices of order 3 (Seberry–Yamada 1990, Theorem 4) gives four ±1 sequences on the group
ℤ₃ × ℤ₁₄₆₁ of order 4383 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The
construction and its proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 6; sign normalisation c = +1.
- T-matrices of order 3 (from base sequences BS(2, 1), built by the program): first rows `+00`, `0+0`, `00+`, `000`.
- Row sums of the four sequences: 107, −67, −37, 15 (107² + (−67)² + (−37)² + 15² = 17532).
- Since gcd(s, t) = 3, the sequences live on ℤ_3 × ℤ_1461, not on ℤ_4383. `sequences.txt` lists the entry at (i, a)
  in position i · 1461 + a, and the Goethals–Seidel array uses the group-developed matrices and the reversal x ↦ −x
  of the group.

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 5843 3` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 4383 has t = 3: the smallest known Hadamard matrix with
odd part 4383 had order 2³ · 4383. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 5844 is equivalent
  to an array of sixteen circulant blocks of order 1461 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and
  Yamada (1990, Theorem 4) multiply such an array by T-matrices of order 3.
- Paley's constructions do not apply: 4n − 1 = 17531 = 47 · 373 and 2n − 1 = 8765 = 5 · 1753 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4); Yamada–Kiyasu does not apply: n − 2 = 4381 = 13 · 337 is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(17532, existence=True)`
  returns Unknown.

We found no construction or statement of order 17532. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 4383 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 17532
```

SHA-256 of `H17532.csv`: `62a113ab4b52153598387c3c0e26349880cdf86820465b5a18e7d7382f26330d`

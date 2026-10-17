# Hadamard matrix of order 23724

23724 = 4 · 5931, with 5931 = 3² · 659.

## Construction

Paley's skew-Hadamard matrix of order 7908 (q = 7907, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 1977 once its points are listed along the four orbits of a torus element of PSL(2, 7907) and rescaled by
signs. Multiplying it by T-matrices of order 3 (Seberry–Yamada 1990, Theorem 4) gives four ±1 sequences on the group
ℤ₃ × ℤ₁₉₇₇ of order 5931 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The
construction and its proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 3 (from base sequences BS(2, 1), built by the program): first rows `+00`, `0+0`, `00+`, `000`.
- Row sums of the four sequences: 1, −99, 101, −61 (1² + (−99)² + 101² + (−61)² = 23724).
- Since gcd(s, t) = 3, the sequences live on ℤ_3 × ℤ_1977, not on ℤ_5931. `sequences.txt` lists the entry at (i, a)
  in position i · 1977 + a, and the Goethals–Seidel array uses the group-developed matrices and the reversal x ↦ −x
  of the group.

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 7907 3` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 5931 has t = 3: the smallest known Hadamard matrix with
odd part 5931 had order 2³ · 5931. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 7908 is equivalent
  to an array of sixteen circulant blocks of order 1977 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and
  Yamada (1990, Theorem 4) multiply such an array by T-matrices of order 3.
- Paley's constructions do not apply: 4n − 1 = 23723 = 7 · 3389 and 2n − 1 = 11861 = 29 · 409 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4); Yamada–Kiyasu does not apply: n − 2 = 5929 = 7² · 11² is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(23724, existence=True)`
  returns Unknown.

We found no construction or statement of order 23724. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 5931 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 23724
```

SHA-256 of `H23724.csv`: `c714a7b7b0a9f42fa8c4f4b04a1735f13941e021ec90d5173be44d7368096455`

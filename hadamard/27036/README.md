# Hadamard matrix of order 27036

27036 = 4 · 6759, with 6759 = 3² · 751.

## Construction

Paley's skew-Hadamard matrix of order 9012 (q = 9011, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 2253 once its points are listed along the four orbits of a torus element of PSL(2, 9011) and rescaled by
signs. Multiplying it by T-matrices of order 3 (Seberry–Yamada 1990, Theorem 4) gives four ±1 sequences on the group
ℤ₃ × ℤ₂₂₅₃ of order 6759 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The
construction and its proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 5; sign normalisation c = +1.
- T-matrices of order 3 (from base sequences BS(2, 1), built by the program): first rows `+00`, `0+0`, `00+`, `000`.
- Row sums of the four sequences: 1, −115, 117, −11 (1² + (−115)² + 117² + (−11)² = 27036).
- Since gcd(s, t) = 3, the sequences live on ℤ_3 × ℤ_2253, not on ℤ_6759. `sequences.txt` lists the entry at (i, a)
  in position i · 2253 + a, and the Goethals–Seidel array uses the group-developed matrices and the reversal x ↦ −x
  of the group.

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 9011 3` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 6759 has t = 8: the smallest known Hadamard matrix with
odd part 6759 had order 2⁸ · 6759. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 9012 is equivalent
  to an array of sixteen circulant blocks of order 2253 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and
  Yamada (1990, Theorem 4) multiply such an array by T-matrices of order 3.
- Paley's constructions do not apply: 4n − 1 = 27035 = 5 · 5407 and 2n − 1 = 13517 = 7 · 1931 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4); Yamada–Kiyasu does not apply: n − 2 = 6757 = 29 · 233 is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(27036, existence=True)`
  returns Unknown.

We found no construction or statement of order 27036. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 6759 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 27036
```

SHA-256 of `H27036.csv`: `6ebfdec4b8d97fcd57aabdb6de9d028e46d86a69ccbec527554621e87b01362d`

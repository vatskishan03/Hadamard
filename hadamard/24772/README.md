# Hadamard matrix of order 24772

24772 = 4 · 6193, with 6193 = 11 · 563.

## Construction

Paley's skew-Hadamard matrix of order 2252 (q = 2251, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 563 once its points are listed along the four orbits of a torus element of PSL(2, 2251) and rescaled by
signs. Multiplying it by T-matrices of order 11 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 6193 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = −1.
- T-matrices of order 11 ([`t-matrices/T11.txt`](../../t-matrices/T11.txt), T-sequences of length 11 from a Golay
  pair of length 10 (Turyn 1974)): first rows `+0000000000`, `0++00+0+-00`, `000++0-00-+`, `00000000000`.
- Row sums of the four sequences: −89, 61, −83, 79 ((−89)² + 61² + (−83)² + 79² = 24772).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 2251 11
t-matrices/T11.txt` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 6193 has t = 7: the smallest known Hadamard matrix with
odd part 6193 had order 2⁷ · 6193. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 2252 is equivalent
  to an array of sixteen circulant blocks of order 563 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 11.
- Paley's constructions do not apply: 4n − 1 = 24771 = 3 · 23 · 359 and 2n − 1 = 12385 = 5 · 2477 are not prime
  powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 6194 cannot exist (6193 is not a sum
  of two squares); Yamada–Kiyasu does not apply: n − 2 = 6191 = 41 · 151 is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(24772, existence=True)`
  returns Unknown.

We found no construction or statement of order 24772. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 6193 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 24772
```

SHA-256 of `H24772.csv`: `b1d79c25163dcc8dc20e9f30439b96e586424f02d2274ad4c1d5240f915b251d`

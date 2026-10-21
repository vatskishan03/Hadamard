# Hadamard matrix of order 32692

32692 = 4 · 8173, with 8173 = 11 · 743.

## Construction

Paley's skew-Hadamard matrix of order 2972 (q = 2971, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 743 once its points are listed along the four orbits of a torus element of PSL(2, 2971) and rescaled by
signs. Multiplying it by T-matrices of order 11 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 8173 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = +1.
- T-matrices of order 11 ([`t-matrices/T11.txt`](../../t-matrices/T11.txt), T-sequences of length 11 from a Golay
  pair of length 10 (Turyn 1974)): first rows `+0000000000`, `0++00+0+-00`, `000++0-00-+`, `00000000000`.
- Row sums of the four sequences: 19, −59, 169, −17 (19² + (−59)² + 169² + (−17)² = 32692).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 2971 11
t-matrices/T11.txt` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 8173 has t = 3: the smallest known Hadamard matrix with
odd part 8173 had order 2³ · 8173. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 2972 is equivalent
  to an array of sixteen circulant blocks of order 743 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 11.
- Paley's constructions do not apply: 4n − 1 = 32691 = 3 · 17 · 641 and 2n − 1 = 16345 = 5 · 7 · 467 are not prime
  powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 8174 cannot exist (8173 is not a sum
  of two squares); Yamada–Kiyasu does not apply: n − 2 = 8171 is ≡ 3 (mod 8).
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(32692, existence=True)`
  returns Unknown.

We found no construction or statement of order 32692. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 8173 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 32692
```

SHA-256 of `H32692.csv`: `5e2ae36315ea14b519b436946a8c394a9e023216403dc69668b75f0c454c25e6`

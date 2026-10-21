# Hadamard matrix of order 30508

30508 = 4 · 7627, with 7627 = 29 · 263.

## Construction

Paley's skew-Hadamard matrix of order 1052 (q = 1051, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 263 once its points are listed along the four orbits of a torus element of PSL(2, 1051) and rescaled by
signs. Multiplying it by T-matrices of order 29 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 7627 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = −1.
- T-matrices of order 29 ([`t-matrices/T29.txt`](../../t-matrices/T29.txt), T-sequences of length 29 as stored in
  SageMath's `t_sequences` module): first rows `++00++-0-++00+000000000000000`, `00-+000+000+-0+00000000000000`,
  `000000000000000+000-0000+000-`, `0000000000000000+++0-+-+0---0`.
- Row sums of the four sequences: 31, −63, 63, 147 (31² + (−63)² + 63² + 147² = 30508).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 1051 29
t-matrices/T29.txt` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 7627 has t = 3: the smallest known Hadamard matrix with
odd part 7627 had order 2³ · 7627. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 1052 is equivalent
  to an array of sixteen circulant blocks of order 263 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 29.
- Paley's constructions do not apply: 4n − 1 = 30507 = 3 · 10169 and 2n − 1 = 15253 = 7 · 2179 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4); Yamada–Kiyasu does not apply: n − 2 = 7625 = 5³ · 61 is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(30508, existence=True)`
  returns Unknown.

We found no construction or statement of order 30508. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 7627 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 30508
```

SHA-256 of `H30508.csv`: `80d400bbb7a516529de340c01384b1ab79224e9d924e4c0dc4d91cf24ed728e9`

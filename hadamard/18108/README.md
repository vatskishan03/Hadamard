# Hadamard matrix of order 18108

18108 = 4 · 4527, with 4527 = 3² · 503.

## Construction

Paley's skew-Hadamard matrix of order 2012 (q = 2011, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 503 once its points are listed along the four orbits of a torus element of PSL(2, 2011) and rescaled by
signs. Multiplying it by T-matrices of order 9 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 4527 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = +1.
- T-matrices of order 9 ([`t-matrices/T9.txt`](../../t-matrices/T9.txt), T-sequences of length 9 from a Golay pair of
  length 8 (Turyn 1974)): first rows `+00000000`, `0+++-0000`, `00000++-+`, `000000000`.
- Row sums of the four sequences: 69, −59, 29, 95 (69² + (−59)² + 29² + 95² = 18108).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 2011 9
t-matrices/T9.txt` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 4527 has t = 3: the smallest known Hadamard matrix with
odd part 4527 had order 2³ · 4527. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 2012 is equivalent
  to an array of sixteen circulant blocks of order 503 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 9.
- Paley's constructions do not apply: 4n − 1 = 18107 = 19 · 953 and 2n − 1 = 9053 = 11 · 823 are not prime powers.
- Miyamoto and Yamada–Kiyasu do not apply: 4527 is not a prime power, and neither is n − 2 = 4525 = 5² · 181.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(18108, existence=True)`
  returns Unknown.

We found no construction or statement of order 18108. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 4527 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 18108
```

SHA-256 of `H18108.csv`: `345697113a58d9b3caf364f51efa7bfd38eabb995211098eb13043eec848795f`

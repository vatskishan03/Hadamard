# Hadamard matrix of order 22132

22132 = 4 · 5533, with 5533 = 11 · 503.

## Construction

Paley's skew-Hadamard matrix of order 2012 (q = 2011, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 503 once its points are listed along the four orbits of a torus element of PSL(2, 2011) and rescaled by
signs. Multiplying it by T-matrices of order 11 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 5533 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = +1.
- T-matrices of order 11 ([`t-matrices/T11.txt`](../../t-matrices/T11.txt), T-sequences of length 11 from a Golay
  pair of length 10 (Turyn 1974)): first rows `+0000000000`, `0++00+0+-00`, `000++0-00-+`, `00000000000`.
- Row sums of the four sequences: 121, −49, 37, 61 (121² + (−49)² + 37² + 61² = 22132).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 2011 11
t-matrices/T11.txt` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 5533 has t = 3: the smallest known Hadamard matrix with
odd part 5533 had order 2³ · 5533. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 2012 is equivalent
  to an array of sixteen circulant blocks of order 503 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 11.
- Paley's constructions do not apply: 4n − 1 = 22131 = 3² · 2459 and 2n − 1 = 11065 = 5 · 2213 are not prime powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 5534 cannot exist (5533 is not a sum
  of two squares); Yamada–Kiyasu does not apply: n − 2 = 5531 is ≡ 3 (mod 8).
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(22132, existence=True)`
  returns Unknown.

We found no construction or statement of order 22132. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 5533 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 22132
```

SHA-256 of `H22132.csv`: `687cc399a05334ab98bebe2334384fe96fe895d39c97ef690950bb1d3f71cf62`

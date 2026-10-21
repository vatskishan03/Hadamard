# Hadamard matrix of order 34204

34204 = 4 · 8551, with 8551 = 17 · 503.

## Construction

Paley's skew-Hadamard matrix of order 2012 (q = 2011, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 503 once its points are listed along the four orbits of a torus element of PSL(2, 2011) and rescaled by
signs. Multiplying it by T-matrices of order 17 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 8551 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = +1.
- T-matrices of order 17 ([`t-matrices/T17.txt`](../../t-matrices/T17.txt), T-sequences of length 17 from a Golay
  pair of length 16 (Turyn 1974)): first rows `+0000000000000000`, `0+++-++-+00000000`, `000000000+++---+-`,
  `00000000000000000`.
- Row sums of the four sequences: 173, −39, 45, 27 (173² + (−39)² + 45² + 27² = 34204).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 2011 17
t-matrices/T17.txt` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 8551 has t = 3: the smallest known Hadamard matrix with
odd part 8551 had order 2³ · 8551. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 2012 is equivalent
  to an array of sixteen circulant blocks of order 503 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 17.
- Paley's constructions do not apply: 4n − 1 = 34203 = 3 · 13 · 877 and 2n − 1 = 17101 = 7² · 349 are not prime
  powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4); Yamada–Kiyasu does not apply: n − 2 = 8549 = 83 · 103 is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(34204, existence=True)`
  returns Unknown.

We found no construction or statement of order 34204. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 8551 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 34204
```

SHA-256 of `H34204.csv`: `dc61b29f2ac565a3f779211154b21b2579c65a2bca4700860a47828e19aa52d9`

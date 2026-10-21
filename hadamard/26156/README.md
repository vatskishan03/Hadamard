# Hadamard matrix of order 26156

26156 = 4 · 6539, with 6539 = 13 · 503.

## Construction

Paley's skew-Hadamard matrix of order 2012 (q = 2011, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 503 once its points are listed along the four orbits of a torus element of PSL(2, 2011) and rescaled by
signs. Multiplying it by T-matrices of order 13 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 6539 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = +1.
- T-matrices of order 13 (from base sequences BS(7, 6), built by the program): first rows `+0+0+00000000`,
  `0+0+0-+000000`, `0000000+0000-`, `00000000++--0`.
- Row sums of the four sequences: 89, −127, 45, −9 (89² + (−127)² + 45² + (−9)² = 26156).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 2011 13` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 6539 has t = 4: the smallest known Hadamard matrix with
odd part 6539 had order 2⁴ · 6539. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 2012 is equivalent
  to an array of sixteen circulant blocks of order 503 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 13.
- Paley's constructions do not apply: 4n − 1 = 26155 = 5 · 5231 and 2n − 1 = 13077 = 3² · 1453 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4); Yamada–Kiyasu does not apply: n − 2 = 6537 = 3 · 2179 is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(26156, existence=True)`
  returns Unknown.

We found no construction or statement of order 26156. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 6539 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 26156
```

SHA-256 of `H26156.csv`: `f3ba5d3eba30bb683bcef550deb09cc4bfa2deda8d439af6b3364cdbeb95ed48`

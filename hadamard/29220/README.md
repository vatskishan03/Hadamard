# Hadamard matrix of order 29220

29220 = 4 · 7305, with 7305 = 3 · 5 · 487.

## Construction

Paley's skew-Hadamard matrix of order 9740 (q = 9739, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 2435 once its points are listed along the four orbits of a torus element of PSL(2, 9739) and rescaled by
signs. Multiplying it by T-matrices of order 3 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 7305 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = +1.
- T-matrices of order 3 (from base sequences BS(2, 1), built by the program): first rows `+00`, `0+0`, `00+`, `000`.
- Row sums of the four sequences: 119, −53, −63, −91 (119² + (−53)² + (−63)² + (−91)² = 29220).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 9739 3` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 7305 has t = 3: the smallest known Hadamard matrix with
odd part 7305 had order 2³ · 7305. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 9740 is equivalent
  to an array of sixteen circulant blocks of order 2435 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and
  Yamada (1990, Theorem 4) multiply such an array by T-matrices of order 3.
- Paley's constructions do not apply: 4n − 1 = 29219 = 61 · 479 and 2n − 1 = 14609 = 7 · 2087 are not prime powers.
- Miyamoto and Yamada–Kiyasu do not apply: 7305 is not a prime power, and neither is n − 2 = 7303 = 67 · 109.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(29220, existence=True)`
  returns Unknown.

We found no construction or statement of order 29220. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 7305 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 29220
```

SHA-256 of `H29220.csv`: `085c2f62b74181c66fef208cccba45bd172306e631bb623ccfd11bbde29d5e18`

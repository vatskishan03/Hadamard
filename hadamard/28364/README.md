# Hadamard matrix of order 28364

28364 = 4 · 7091, with 7091 = 7 · 1013.

## Construction

Paley's skew-Hadamard matrix of order 4052 (q = 4051, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 1013 once its points are listed along the four orbits of a torus element of PSL(2, 4051) and rescaled by
signs. Multiplying it by T-matrices of order 7 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 7091 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = +1.
- T-matrices of order 7 (from base sequences BS(4, 3), built by the program): first rows `+000000`, `0++-000`,
  `0000+0+`, `00000+0`.
- Row sums of the four sequences: 3, −43, 91, −135 (3² + (−43)² + 91² + (−135)² = 28364).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 4051 7` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 7091 has t = 3: the smallest known Hadamard matrix with
odd part 7091 had order 2³ · 7091. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 4052 is equivalent
  to an array of sixteen circulant blocks of order 1013 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and
  Yamada (1990, Theorem 4) multiply such an array by T-matrices of order 7.
- Paley's constructions do not apply: 4n − 1 = 28363 = 113 · 251 and 2n − 1 = 14181 = 3 · 29 · 163 are not prime
  powers.
- Miyamoto and Yamada–Kiyasu do not apply: 7091 is not a prime power, and neither is n − 2 = 7089 = 3 · 17 · 139.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(28364, existence=True)`
  returns Unknown.

We found no construction or statement of order 28364. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 7091 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 28364
```

SHA-256 of `H28364.csv`: `3da6576b1f1ea8cfb752d49d06770df392e40132e78f023417c431106d692014`

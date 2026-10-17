# Hadamard matrix of order 36596

36596 = 4 · 9149, with 9149 = 7 · 1307.

## Construction

Paley's skew-Hadamard matrix of order 5228 (q = 5227, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 1307 once its points are listed along the four orbits of a torus element of PSL(2, 5227) and rescaled by
signs. Multiplying it by T-matrices of order 7 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 9149 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 7 (from base sequences BS(4, 3), built by the program): first rows `+000000`, `0++-000`,
  `0000+0+`, `00000+0`.
- Row sums of the four sequences: −91, −55, −3, 159 ((−91)² + (−55)² + (−3)² + 159² = 36596).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 5227 7` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 9149 has t = 3: the smallest known Hadamard matrix with
odd part 9149 had order 2³ · 9149. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 5228 is equivalent
  to an array of sixteen circulant blocks of order 1307 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and
  Yamada (1990, Theorem 4) multiply such an array by T-matrices of order 7.
- Paley's constructions do not apply: 4n − 1 = 36595 = 5 · 13 · 563 and 2n − 1 = 18297 = 3² · 19 · 107 are not prime
  powers.
- Miyamoto and Yamada–Kiyasu do not apply: 9149 is not a prime power, and neither is n − 2 = 9147 = 3 · 3049.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(36596, existence=True)`
  returns Unknown.

We found no construction or statement of order 36596. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 9149 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 36596
```

SHA-256 of `H36596.csv`: `844c1fbba26511d8dca878783a40601deb676bec0fb14e1ea63c36c47b31b427`

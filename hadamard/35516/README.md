# Hadamard matrix of order 35516

35516 = 4 · 8879, with 8879 = 13 · 683.

## Construction

Paley's skew-Hadamard matrix of order 2732 (q = 2731, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 683 once its points are listed along the four orbits of a torus element of PSL(2, 2731) and rescaled by
signs. Multiplying it by T-matrices of order 13 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 8879 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = +1.
- T-matrices of order 13 (from base sequences BS(7, 6), built by the program): first rows `+0+0+00000000`,
  `0+0+0-+000000`, `0000000+0000-`, `00000000++--0`.
- Row sums of the four sequences: 45, −61, −149, 87 (45² + (−61)² + (−149)² + 87² = 35516).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 2731 13` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 8879 has t = 3: the smallest known Hadamard matrix with
odd part 8879 had order 2³ · 8879. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 2732 is equivalent
  to an array of sixteen circulant blocks of order 683 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 13.
- Paley's constructions do not apply: 4n − 1 = 35515 = 5 · 7103 and 2n − 1 = 17757 = 3² · 1973 are not prime powers.
- Miyamoto and Yamada–Kiyasu do not apply: 8879 is not a prime power, and neither is n − 2 = 8877 = 3 · 11 · 269.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(35516, existence=True)`
  returns Unknown.

We found no construction or statement of order 35516. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 8879 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 35516
```

SHA-256 of `H35516.csv`: `f9665e6b360c75847513b6082c4fd73686bd97569c9409ddecc2082b4c996cdc`

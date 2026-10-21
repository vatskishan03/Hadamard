# Hadamard matrix of order 36204

36204 = 4 · 9051, with 9051 = 3 · 7 · 431.

## Construction

Paley's skew-Hadamard matrix of order 5172 (q = 5171, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 1293 once its points are listed along the four orbits of a torus element of PSL(2, 5171) and rescaled by
signs. Multiplying it by T-matrices of order 7 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 9051 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 6; sign normalisation c = −1.
- T-matrices of order 7 (from base sequences BS(4, 3), built by the program): first rows `+000000`, `0++-000`,
  `0000+0+`, `00000+0`.
- Row sums of the four sequences: −51, 89, 59, −149 ((−51)² + 89² + 59² + (−149)² = 36204).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 5171 7` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 9051 has t = 3: the smallest known Hadamard matrix with
odd part 9051 had order 2³ · 9051. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 5172 is equivalent
  to an array of sixteen circulant blocks of order 1293 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and
  Yamada (1990, Theorem 4) multiply such an array by T-matrices of order 7.
- Paley's constructions do not apply: 4n − 1 = 36203 = 41 · 883 and 2n − 1 = 18101 = 23 · 787 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4). Yamada's theorem with q = n − 2 = 9049 needs a symmetric conference matrix of order 4526; 4525 = 5² · 181 is
  not a prime power and we know no such matrix.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(36204, existence=True)`
  returns Unknown.

We found no construction or statement of order 36204. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 9051 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 36204
```

SHA-256 of `H36204.csv`: `d756b8d73500de767122e31d0b038e51e91d489a3397ce6de4c438ee203ca47c`

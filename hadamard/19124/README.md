# Hadamard matrix of order 19124

19124 = 4 · 4781, with 4781 = 7 · 683.

## Construction

Paley's skew-Hadamard matrix of order 2732 (q = 2731, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 683 once its points are listed along the four orbits of a torus element of PSL(2, 2731) and rescaled by
signs. Multiplying it by T-matrices of order 7 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 4781 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = +1.
- T-matrices of order 7 (from base sequences BS(4, 3), built by the program): first rows `+000000`, `0++-000`,
  `0000+0+`, `00000+0`.
- Row sums of the four sequences: 21, 87, −83, 65 (21² + 87² + (−83)² + 65² = 19124).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 2731 7` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 4781 has t = 3: the smallest known Hadamard matrix with
odd part 4781 had order 2³ · 4781. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 2732 is equivalent
  to an array of sixteen circulant blocks of order 683 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 7.
- Paley's constructions do not apply: 4n − 1 = 19123 = 13 · 1471 and 2n − 1 = 9561 = 3 · 3187 are not prime powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 4782 cannot exist (4781 is not a sum
  of two squares); Yamada–Kiyasu does not apply: n − 2 = 4779 = 3⁴ · 59 is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(19124, existence=True)`
  returns Unknown.

We found no construction or statement of order 19124. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 4781 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 19124
```

SHA-256 of `H19124.csv`: `fb97b9ae483caf0e40319949a492f30b781c47ac855151570247803b5a35515f`

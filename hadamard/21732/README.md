# Hadamard matrix of order 21732

21732 = 4 · 5433, with 5433 = 3 · 1811.

## Construction

Paley's skew-Hadamard matrix of order 7244 (q = 7243, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 1811 once its points are listed along the four orbits of a torus element of PSL(2, 7243) and rescaled by
signs. Multiplying it by T-matrices of order 3 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 5433 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 3 (from base sequences BS(2, 1), built by the program): first rows `+00`, `0+0`, `00+`, `000`.
- Row sums of the four sequences: −87, 1, 89, 79 ((−87)² + 1² + 89² + 79² = 21732).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 7243 3` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 5433 has t = 4: the smallest known Hadamard matrix with
odd part 5433 had order 2⁴ · 5433. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 7244 is equivalent
  to an array of sixteen circulant blocks of order 1811 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and
  Yamada (1990, Theorem 4) multiply such an array by T-matrices of order 3.
- Paley's constructions do not apply: 4n − 1 = 21731 = 31 · 701 and 2n − 1 = 10865 = 5 · 41 · 53 are not prime
  powers.
- Miyamoto and Yamada–Kiyasu do not apply: 5433 is not a prime power, and n − 2 = 5431 is ≡ 7 (mod 8).
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(21732, existence=True)`
  returns Unknown.

We found no construction or statement of order 21732. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 5433 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 21732
```

SHA-256 of `H21732.csv`: `c66f8b921552d08c3e87ae759dfd25afca7555019a98c7c72986899139a1e623`

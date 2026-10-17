# Hadamard matrix of order 13092

13092 = 4 · 3273, with 3273 = 3 · 1091.

## Construction

Paley's skew-Hadamard matrix of order 4364 (q = 4363, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 1091 once its points are listed along the four orbits of a torus element of PSL(2, 4363) and rescaled by
signs. Multiplying it by T-matrices of order 3 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 3273 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 3 (from base sequences BS(2, 1), built by the program): first rows `+00`, `0+0`, `00+`, `000`.
- Row sums of the four sequences: −37, −41, 81, −59 ((−37)² + (−41)² + 81² + (−59)² = 13092).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 4363 3` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 3273 has t = 3: the smallest known Hadamard matrix with
odd part 3273 had order 2³ · 3273. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 4364 is equivalent
  to an array of sixteen circulant blocks of order 1091 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and
  Yamada (1990, Theorem 4) multiply such an array by T-matrices of order 3.
- Paley's constructions do not apply: 4n − 1 = 13091 = 13 · 19 · 53 and 2n − 1 = 6545 = 5 · 7 · 11 · 17 are not prime
  powers.
- Miyamoto and Yamada–Kiyasu do not apply: 3273 is not a prime power, and n − 2 = 3271 is ≡ 7 (mod 8).
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(13092, existence=True)`
  returns Unknown.

We found no construction or statement of order 13092. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 3273 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 13092
```

SHA-256 of `H13092.csv`: `e36f8a0f00bd260ff35083c5b87d4f775382833d567c0995d37e4daf6941ab29`

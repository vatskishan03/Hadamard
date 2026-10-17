# Hadamard matrix of order 39012

39012 = 4 · 9753, with 9753 = 3 · 3251.

## Construction

Paley's skew-Hadamard matrix of order 13004 (q = 13003, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 3251 once its points are listed along the four orbits of a torus element of PSL(2, 13003) and rescaled by
signs. Multiplying it by T-matrices of order 3 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 9753 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 3 (from base sequences BS(2, 1), built by the program): first rows `+00`, `0+0`, `00+`, `000`.
- Row sums of the four sequences: 139, −127, −9, −59 (139² + (−127)² + (−9)² + (−59)² = 39012).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 13003 3` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 9753 has t = 3: the smallest known Hadamard matrix with
odd part 9753 had order 2³ · 9753. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 13004 is equivalent
  to an array of sixteen circulant blocks of order 3251 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and
  Yamada (1990, Theorem 4) multiply such an array by T-matrices of order 3.
- Paley's constructions do not apply: 4n − 1 = 39011 = 7 · 5573 and 2n − 1 = 19505 = 5 · 47 · 83 are not prime
  powers.
- Miyamoto and Yamada–Kiyasu do not apply: 9753 is not a prime power, and neither is n − 2 = 9751 = 7² · 199.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(39012, existence=True)`
  returns Unknown.

We found no construction or statement of order 39012. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 9753 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 39012
```

SHA-256 of `H39012.csv`: `486f118c2188e9fa229a6ff0c2b4d5ef46a5e321dcfb2688abdc5bed2d8ab4cf`

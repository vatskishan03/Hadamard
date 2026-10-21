# Hadamard matrix of order 12052

12052 = 4 · 3013, with 3013 = 23 · 131.

## Construction

Paley's skew-Hadamard matrix of order 524 (q = 523, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of order
s = 131 once its points are listed along the four orbits of a torus element of PSL(2, 523) and rescaled by signs.
Multiplying it by T-matrices of order 23 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of length
3013 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its proof
are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 23 ([`t-matrices/T23.txt`](../../t-matrices/T23.txt), T-sequences of length 23 as stored in
  SageMath's `t_sequences` module): first rows `+0000000000000000000000`, `0+++++--+-+-00000000000`,
  `000000000000+0+0-0+0-0+`, `0000000000000+0+0-0+0+0`.
- Row sums of the four sequences: 17, −81, 21, 69 (17² + (−81)² + 21² + 69² = 12052).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 523 23
t-matrices/T23.txt` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 3013 has t = 3: the smallest known Hadamard matrix with
odd part 3013 had order 2³ · 3013. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 524 is equivalent to
  an array of sixteen circulant blocks of order 131 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 23.
- Paley's constructions do not apply: 4n − 1 = 12051 = 3² · 13 · 103 and 2n − 1 = 6025 = 5² · 241 are not prime
  powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 3014 cannot exist (3013 is not a sum
  of two squares); Yamada–Kiyasu does not apply: n − 2 = 3011 is ≡ 3 (mod 8).
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(12052, existence=True)`
  returns Unknown.

A second route reaches the same order: Cooper–Wallis with T-matrices of order 131 and Williamson-type matrices of
order 23 (Seberry–Yamada 1992, Table A.1). T-sequences of length 131 follow from the Turyn-type sequences TT(44) of
S. London and I. Kotsireas (Cryptogr. Commun. 17 (2025) 1601–1610), which we could not read; if that paper lists
Hadamard orders, this one may be among them. Our T-matrices of order 131 are in
[`t-matrices/T131.txt`](../../t-matrices/T131.txt).

We found no construction or statement of order 12052. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 3013 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 12052
```

SHA-256 of `H12052.csv`: `54c723bda408b0fc99aec77c20bfea27f09070217bb0f8a968dcce73d4c6903d`

# Hadamard matrix of order 28820

28820 = 4 · 7205, with 7205 = 5 · 11 · 131.

## Construction

Paley's skew-Hadamard matrix of order 524 (q = 523, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of order
s = 131 once its points are listed along the four orbits of a torus element of PSL(2, 523) and rescaled by signs.
Multiplying it by T-matrices of order 55 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of length
7205 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its proof
are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 55 ([`t-matrices/T55.txt`](../../t-matrices/T55.txt), from base sequences BS(6, 5) by Yang's
  multiplication by 5 (Seberry–Yamada 1992, Theorem 5.23, with the Golay pair of length 2; `search/yang5.py 5
  --turyn`)): first rows in `certificate.json`.
- Row sums of the four sequences: −167, −9, 25, −15 ((−167)² + (−9)² + 25² + (−15)² = 28820).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 523 55
t-matrices/T55.txt` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 7205 has t = 3: the smallest known Hadamard matrix with
odd part 7205 had order 2³ · 7205. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 524 is equivalent to
  an array of sixteen circulant blocks of order 131 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 55.
- Paley's constructions do not apply: 4n − 1 = 28819 = 7 · 23 · 179 and 2n − 1 = 14409 = 3² · 1601 are not prime
  powers.
- Miyamoto and Yamada–Kiyasu do not apply: 7205 is not a prime power, and neither is n − 2 = 7203 = 3 · 7⁴.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(28820, existence=True)`
  returns Unknown.

A second route reaches the same order: Cooper–Wallis with T-matrices of order 131 and Williamson-type matrices of
order 55 (Seberry–Yamada 1992, Table A.1). T-sequences of length 131 follow from the Turyn-type sequences TT(44) of
S. London and I. Kotsireas (Cryptogr. Commun. 17 (2025) 1601–1610), which we could not read; if that paper lists
Hadamard orders, this one may be among them. Our T-matrices of order 131 are in
[`t-matrices/T131.txt`](../../t-matrices/T131.txt).

We found no construction or statement of order 28820. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 7205 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 28820
```

SHA-256 of `H28820.csv`: `bb4448668a18afbc6f95c4688f5a176da757730a53485dc04fc15bed72a77924`

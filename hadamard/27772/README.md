# Hadamard matrix of order 27772

27772 = 4 · 6943, with 6943 = 53 · 131.

## Construction

Paley's skew-Hadamard matrix of order 524 (q = 523, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of order
s = 131 once its points are listed along the four orbits of a torus element of PSL(2, 523) and rescaled by signs.
Multiplying it by T-matrices of order 53 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of length
6943 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its proof
are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 53 ([`t-matrices/T53.txt`](../../t-matrices/T53.txt), T-sequences of length 53 as stored in
  SageMath's `t_sequences` module): first rows in `certificate.json`.
- Row sums of the four sequences: −101, −89, 91, 37 ((−101)² + (−89)² + 91² + 37² = 27772).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 523 53
t-matrices/T53.txt` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 6943 has t = 8: the smallest known Hadamard matrix with
odd part 6943 had order 2⁸ · 6943. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 524 is equivalent to
  an array of sixteen circulant blocks of order 131 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 53.
- Paley's constructions do not apply: 4n − 1 = 27771 = 3 · 9257 and 2n − 1 = 13885 = 5 · 2777 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4); Yamada–Kiyasu does not apply: n − 2 = 6941 = 11 · 631 is not a prime power.
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(27772, existence=True)`
  returns Unknown.

A second route reaches the same order: Cooper–Wallis with T-matrices of order 131 and Williamson-type matrices of
order 53 (Seberry–Yamada 1992, Table A.1). T-sequences of length 131 follow from the Turyn-type sequences TT(44) of
S. London and I. Kotsireas (Cryptogr. Commun. 17 (2025) 1601–1610), which we could not read; if that paper lists
Hadamard orders, this one may be among them. Our T-matrices of order 131 are in
[`t-matrices/T131.txt`](../../t-matrices/T131.txt).

We found no construction or statement of order 27772. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 6943 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 27772
```

SHA-256 of `H27772.csv`: `526d8e6780d6342b9bec4800d1afb19942cc5b7aa904f0ab0d976086a08e76da`

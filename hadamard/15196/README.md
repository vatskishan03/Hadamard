# Hadamard matrix of order 15196

15196 = 4 · 3799, with 3799 = 29 · 131.

## Construction

Paley's skew-Hadamard matrix of order 524 (q = 523, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of order
s = 131 once its points are listed along the four orbits of a torus element of PSL(2, 523) and rescaled by signs.
Multiplying it by T-matrices of order 29 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of length
3799 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its proof
are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 29 ([`t-matrices/T29.txt`](../../t-matrices/T29.txt), T-sequences of length 29 as stored in
  SageMath's `t_sequences` module): first rows `++00++-0-++00+000000000000000`, `00-+000+000+-0+00000000000000`,
  `000000000000000+000-0000+000-`, `0000000000000000+++0-+-+0---0`.
- Row sums of the four sequences: 23, −43, 113, −7 (23² + (−43)² + 113² + (−7)² = 15196).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 523 29
t-matrices/T29.txt` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 3799 has t = 3: the smallest known Hadamard matrix with
odd part 3799 had order 2³ · 3799. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 524 is equivalent to
  an array of sixteen circulant blocks of order 131 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 29.
- Paley's constructions do not apply: 4n − 1 = 15195 = 3 · 5 · 1013 and 2n − 1 = 7597 = 71 · 107 are not prime
  powers.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(15196, existence=True)`
  returns Unknown.

A second route reaches the same order: Cooper–Wallis with T-matrices of order 131 and Williamson-type matrices of
order 29 (Seberry–Yamada 1992, Table A.1). T-sequences of length 131 follow from the Turyn-type sequences TT(44) of
S. London and I. Kotsireas (Cryptogr. Commun. 17 (2025) 1601–1610), which we could not read; if that paper lists
Hadamard orders, this one may be among them. Our T-matrices of order 131 are in
[`t-matrices/T131.txt`](../../t-matrices/T131.txt).

Kiyasu's theorem with q = n − 2 = 3797 ≡ 5 (mod 8) reaches it from a skew-Hadamard matrix of order 1900, which this
repository constructs ([skew 1900](../../skew/1900)); the census of May 2025 lists v = 475 as undecided.

We found no construction or statement of order 15196. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 3799 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 15196
```

SHA-256 of `H15196.csv`: `61faacb0d373d8720e7d275454df21f61be04616cfe1156f1367409b7a5003cf`

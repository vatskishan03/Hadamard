# Hadamard matrix of order 20436

20436 = 4 · 5109, with 5109 = 3 · 13 · 131.

## Construction

Paley's skew-Hadamard matrix of order 1572 (q = 1571, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 393 once its points are listed along the four orbits of a torus element of PSL(2, 1571) and rescaled by
signs. Multiplying it by T-matrices of order 13 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 5109 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 6; sign normalisation c = −1.
- T-matrices of order 13 (from base sequences BS(7, 6), built by the program): first rows `+0+0+00000000`,
  `0+0+0-+000000`, `0000000+0000-`, `00000000++--0`.
- Row sums of the four sequences: 13, −13, 127, −63 (13² + (−13)² + 127² + (−63)² = 20436).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 1571 13` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 5109 has t = 3: the smallest known Hadamard matrix with
odd part 5109 had order 2³ · 5109. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 1572 is equivalent
  to an array of sixteen circulant blocks of order 393 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 13.
- Paley's constructions do not apply: 4n − 1 = 20435 = 5 · 61 · 67 and 2n − 1 = 10217 = 17 · 601 are not prime
  powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 5110 cannot exist (5109 is not a sum
  of two squares); Yamada–Kiyasu does not apply: n − 2 = 5107 is ≡ 3 (mod 8).
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(20436, existence=True)`
  returns Unknown.

A second route reaches the same order: Cooper–Wallis with T-matrices of order 131 and Williamson-type matrices of
order 39 (Seberry–Yamada 1992, Table A.1). T-sequences of length 131 follow from the Turyn-type sequences TT(44) of
S. London and I. Kotsireas (Cryptogr. Commun. 17 (2025) 1601–1610), which we could not read; if that paper lists
Hadamard orders, this one may be among them. Our T-matrices of order 131 are in
[`t-matrices/T131.txt`](../../t-matrices/T131.txt).

We found no construction or statement of order 20436. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 5109 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 20436
```

SHA-256 of `H20436.csv`: `7386e0375499c976b59c2c6d267688d249ce0976d5ded1638baee6ded2e981e4`

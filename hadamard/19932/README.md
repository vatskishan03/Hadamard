# Hadamard matrix of order 19932

19932 = 4 · 4983, with 4983 = 3 · 11 · 151.

## Construction

Paley's skew-Hadamard matrix of order 1812 (q = 1811, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 453 once its points are listed along the four orbits of a torus element of PSL(2, 1811) and rescaled by
signs. Multiplying it by T-matrices of order 11 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 4983 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 10; sign normalisation c = +1.
- T-matrices of order 11 ([`t-matrices/T11.txt`](../../t-matrices/T11.txt), T-sequences of length 11 from a Golay
  pair of length 10 (Turyn 1974)): first rows `+0000000000`, `0++00+0+-00`, `000++0-00-+`, `00000000000`.
- Row sums of the four sequences: 81, −49, 77, −71 (81² + (−49)² + 77² + (−71)² = 19932).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 1811 11
t-matrices/T11.txt` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 4983 has t = 3: the smallest known Hadamard matrix with
odd part 4983 had order 2³ · 4983. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 1812 is equivalent
  to an array of sixteen circulant blocks of order 453 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 11.
- Paley's constructions do not apply: 4n − 1 = 19931 = 19 · 1049 and 2n − 1 = 9965 = 5 · 1993 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4); Yamada–Kiyasu does not apply: n − 2 = 4981 = 17 · 293 is not a prime power.
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(19932, existence=True)`
  returns Unknown.

Since 17 September 2026 the order also follows by Cooper–Wallis from T-matrices of order 151 (Schumer, Zenodo,
doi:10.5281/zenodo.22820997) and Williamson-type matrices of order 33 (Seberry–Yamada 1992, Table A.1), but we found
no construction or statement of it.

We found no construction or statement of order 19932. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 4983 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 19932
```

SHA-256 of `H19932.csv`: `4f6c4066efa19c0587a641925ac1b9208bdef9bfb6bf9ebe99b2bb41567e0c4e`

# Hadamard matrix of order 32012

32012 = 4 · 8003, with 8003 = 53 · 151.

## Construction

Paley's skew-Hadamard matrix of order 212 (q = 211, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of order
s = 53 once its points are listed along the four orbits of a torus element of PSL(2, 211) and rescaled by signs.
Multiplying it by T-matrices of order 151 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 8003 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = −1.
- T-matrices of order 151 ([`t-matrices/T151.txt`](../../t-matrices/T151.txt), see
  [t-matrices/README.md](../../t-matrices/README.md)): first rows in `certificate.json`.
- Row sums of the four sequences: −111, 29, −9, −137 ((−111)² + 29² + (−9)² + (−137)² = 32012).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 211 151
t-matrices/T151.txt` writes the certificate.

## Status

**New, beyond the tables.** In Table 1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 8003 has t =
3: the smallest known Hadamard matrix with odd part 8003 had order 2³ · 8003. The order lies above the range of
Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 212 is equivalent to
  an array of sixteen circulant blocks of order 53 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 151, which are known only since 2026.
- Paley's constructions do not apply: 4n − 1 = 32011 = 7 · 17 · 269 and 2n − 1 = 16005 = 3 · 5 · 11 · 97 are not
  prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4); Yamada–Kiyasu does not apply: n − 2 = 8001 = 3² · 7 · 127 is not a prime power.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(32012, existence=True)`
  returns Unknown.

Since September 17, 2026 the order also follows by Cooper–Wallis from T-matrices of order 151 (Schumer, Zenodo,
doi:10.5281/zenodo.22820997) and Williamson-type matrices of order 53 (Seberry–Yamada 1992, Table A.1), but we found
no construction or statement of it.

We found no construction or statement of order 32012. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 8003 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 32012
```

SHA-256 of `H32012.csv`: `e6de4af3b548e2343e3fe5ab813ef32a8a1c79ba10f9b959e08e8a61cb74d3c3`

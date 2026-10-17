# Hadamard matrix of order 36844

36844 = 4 · 9211, with 9211 = 61 · 151.

## Construction

Paley's skew-Hadamard matrix of order 244 (q = 243 = 3⁵, a prime power ≡ 3 mod 8) is a 4 × 4 array of circulant
blocks of order s = 61 once its points are listed along the four orbits of a torus element of PSL(2, 243) and
rescaled by signs. Multiplying it by T-matrices of order 151 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1
sequences of length 9211 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The
construction and its proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = x (integer code 3) in GF(243) = GF(3)[x]/(x⁵ + 2x⁴ + 1),
  elements coded as Σ dᵢ 3ⁱ for Σ dᵢ xⁱ; sign normalisation c = +1.
- T-matrices of order 151 ([`t-matrices/T151.txt`](../../t-matrices/T151.txt), see
  [t-matrices/README.md](../../t-matrices/README.md)): first rows in `certificate.json`.
- Row sums of the four sequences: −99, 9, 121, −111 ((−99)² + 9² + 121² + (−111)² = 36844).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 243 151
t-matrices/T151.txt` writes the certificate.

## Status

**New, beyond the tables.** In Table 1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 9211 has t =
3: the smallest known Hadamard matrix with odd part 9211 had order 2³ · 9211. The order lies above the range of
Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 244 is equivalent to
  an array of sixteen circulant blocks of order 61 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 151, which are known only since 2026.
- Paley's constructions do not apply: 4n − 1 = 36843 = 3 · 12281 and 2n − 1 = 18421 = 13² · 109 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4). Yamada's theorem with q = n − 2 = 9209 needs a symmetric conference matrix of order 4606; 4605 = 3 · 5 · 307
  is not a sum of two squares, so no conference matrix of order 4606 exists.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(36844, existence=True)`
  returns Unknown.

Since September 17, 2026 the order also follows by Cooper–Wallis from T-matrices of order 151 (Schumer, Zenodo,
doi:10.5281/zenodo.22820997) and Williamson-type matrices of order 61 (Seberry–Yamada 1992, Table A.1), but we found
no construction or statement of it.

We found no construction or statement of order 36844. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 9211 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 36844
```

SHA-256 of `H36844.csv`: `dc9e3de4103c99cc8417170ba86aa2f996ef62c1698da8e873c2c5f973c03d32`

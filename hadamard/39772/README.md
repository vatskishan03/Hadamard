# Hadamard matrix of order 39772

39772 = 4 · 9943, with 9943 = 61 · 163.

## Construction

Paley's skew-Hadamard matrix of order 244 (q = 243 = 3⁵, a prime power ≡ 3 mod 8) is a 4 × 4 array of circulant
blocks of order s = 61 once its points are listed along the four orbits of a torus element of PSL(2, 243) and
rescaled by signs. Multiplying it by T-matrices of order 163 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1
sequences of length 9943 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The
construction and its proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = x (integer code 3) in GF(243) = GF(3)[x]/(x⁵ + 2x⁴ + 1),
  elements coded as Σ dᵢ 3ⁱ for Σ dᵢ xⁱ; sign normalisation c = +1.
- T-matrices of order 163 ([`t-matrices/T163.txt`](../../t-matrices/T163.txt), see
  [t-matrices/README.md](../../t-matrices/README.md)): first rows in `certificate.json`.
- Row sums of the four sequences: −97, −91, 109, −101 ((−97)² + (−91)² + 109² + (−101)² = 39772).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 243 163
t-matrices/T163.txt` writes the certificate.

## Status

**New, beyond the tables.** In Table 1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 9943 has t =
3: the smallest known Hadamard matrix with odd part 9943 had order 2³ · 9943. The order lies above the range of
Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 244 is equivalent to
  an array of sixteen circulant blocks of order 61 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 163, which are known only since 2026.
- Paley's constructions do not apply: 4n − 1 = 39771 = 3⁴ · 491 and 2n − 1 = 19885 = 5 · 41 · 97 are not prime
  powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4). Kiyasu's theorem with q = n − 2 = 9941 needs a skew-Hadamard matrix of order 4972 = 4 · 1243, and v = 1243 is
  undecided in the Đoković–Seberry census.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(39772, existence=True)`
  returns Unknown.

Since August 29, 2026 the order also follows by Cooper–Wallis from T-matrices of order 163 (JD-Jones-ASES/Hadamard-T)
and Williamson-type matrices of order 61 (Seberry–Yamada 1992, Table A.1), but we found no construction or statement
of it.

We found no construction or statement of order 39772. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 9943 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 39772
```

SHA-256 of `H39772.csv`: `ff6ae86d2db54605391e9470608a479c8ca4984ec503caee9d23086642359785`

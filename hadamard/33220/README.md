# Hadamard matrix of order 33220

33220 = 4 · 8305, with 8305 = 5 · 11 · 151.

## Construction

Paley's skew-Hadamard matrix of order 3020 (q = 3019, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 755 once its points are listed along the four orbits of a torus element of PSL(2, 3019) and rescaled by
signs. Multiplying it by T-matrices of order 11 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 8305 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = +1.
- T-matrices of order 11 ([`t-matrices/T11.txt`](../../t-matrices/T11.txt), T-sequences of length 11 from a Golay
  pair of length 10 (Turyn 1974)): first rows `+0000000000`, `0++00+0+-00`, `000++0-00-+`, `00000000000`.
- Row sums of the four sequences: 55, −59, 133, −95 (55² + (−59)² + 133² + (−95)² = 33220).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 3019 11
t-matrices/T11.txt` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 8305 has t = 3: the smallest known Hadamard matrix with
odd part 8305 had order 2³ · 8305. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 3020 is equivalent
  to an array of sixteen circulant blocks of order 755 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 11.
- Paley's constructions do not apply: 4n − 1 = 33219 = 3² · 3691 and 2n − 1 = 16609 = 17 · 977 are not prime powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 8306 cannot exist (8305 is not a sum
  of two squares); Yamada–Kiyasu does not apply: n − 2 = 8303 = 19² · 23 is not a prime power.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(33220, existence=True)`
  returns Unknown.

Since September 17, 2026 the order also follows by Cooper–Wallis from T-matrices of order 151 (Schumer, Zenodo,
doi:10.5281/zenodo.22820997) and Williamson-type matrices of order 55 (Seberry–Yamada 1992, Table A.1), but we found
no construction or statement of it.

We found no construction or statement of order 33220. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 8305 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 33220
```

SHA-256 of `H33220.csv`: `5548e3accfb85d3e68cdad311ba7f97164f5508a9d25b0911b36ff21d927ecb7`

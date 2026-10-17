# Hadamard matrix of order 34556

34556 = 4 · 8639, with 8639 = 53 · 163.

## Construction

Paley's skew-Hadamard matrix of order 212 (q = 211, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of order
s = 53 once its points are listed along the four orbits of a torus element of PSL(2, 211) and rescaled by signs.
Multiplying it by T-matrices of order 163 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 8639 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = −1.
- T-matrices of order 163 ([`t-matrices/T163.txt`](../../t-matrices/T163.txt), see
  [t-matrices/README.md](../../t-matrices/README.md)): first rows in `certificate.json`.
- Row sums of the four sequences: −63, −45, 1, −169 ((−63)² + (−45)² + 1² + (−169)² = 34556).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 211 163
t-matrices/T163.txt` writes the certificate.

## Status

**New, beyond the tables.** In Table 1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 8639 has t =
3: the smallest known Hadamard matrix with odd part 8639 had order 2³ · 8639. The order lies above the range of
Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 212 is equivalent to
  an array of sixteen circulant blocks of order 53 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 163, which are known only since 2026.
- Paley's constructions do not apply: 4n − 1 = 34555 = 5 · 6911 and 2n − 1 = 17277 = 3 · 13 · 443 are not prime
  powers.
- Miyamoto and Yamada–Kiyasu do not apply: 8639 is not a prime power, and neither is n − 2 = 8637 = 3 · 2879.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(34556, existence=True)`
  returns Unknown.

Since August 29, 2026 the order also follows by Cooper–Wallis from T-matrices of order 163 (JD-Jones-ASES/Hadamard-T)
and Williamson-type matrices of order 53 (Seberry–Yamada 1992, Table A.1), but we found no construction or statement
of it.

We found no construction or statement of order 34556. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 8639 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 34556
```

SHA-256 of `H34556.csv`: `a345b368ba85e46f7014192b750f45d465de8f69b98e324e5c0a2ab2f3da950e`

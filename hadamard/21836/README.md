# Hadamard matrix of order 21836

21836 = 4 · 5459, with 5459 = 53 · 103.

## Construction

Paley's skew-Hadamard matrix of order 212 (q = 211, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of order
s = 53 once its points are listed along the four orbits of a torus element of PSL(2, 211) and rescaled by signs.
Multiplying it by T-matrices of order 103 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 5459 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = −1.
- T-matrices of order 103 ([`t-matrices/T103.txt`](../../t-matrices/T103.txt), see
  [t-matrices/README.md](../../t-matrices/README.md)): first rows in `certificate.json`.
- Row sums of the four sequences: −63, 51, 79, −95 ((−63)² + 51² + 79² + (−95)² = 21836).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 211 103
t-matrices/T103.txt` writes the certificate.

## Status

**New, beyond the tables.** In Table 1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 5459 has t =
3: the smallest known Hadamard matrix with odd part 5459 had order 2³ · 5459. The order lies above the range of
Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 212 is equivalent to
  an array of sixteen circulant blocks of order 53 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 103, which are known only since 2026.
- Paley's constructions do not apply: 4n − 1 = 21835 = 5 · 11 · 397 and 2n − 1 = 10917 = 3² · 1213 are not prime
  powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4); Yamada–Kiyasu does not apply: n − 2 = 5457 = 3 · 17 · 107 is not a prime power.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(21836, existence=True)`
  returns Unknown.

Since August 29, 2026 the order also follows by Cooper–Wallis from T-matrices of order 103 (JD-Jones-ASES/Hadamard-T,
decoded from Schneider's matrix of order 2060) and Williamson-type matrices of order 53 (Seberry–Yamada 1992, Table
A.1), but we found no construction or statement of it.

We found no construction or statement of order 21836. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 5459 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 21836
```

SHA-256 of `H21836.csv`: `92106bbb87a623c733cedd3978b744fb7455f973fbdd0bff1e20e1fd461bbb18`

# Hadamard matrix of order 16892

16892 = 4 · 4223, with 4223 = 41 · 103.

## Construction

Paley's skew-Hadamard matrix of order 164 (q = 163, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of order
s = 41 once its points are listed along the four orbits of a torus element of PSL(2, 163) and rescaled by signs.
Multiplying it by T-matrices of order 103 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 4223 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 103 ([`t-matrices/T103.txt`](../../t-matrices/T103.txt), see
  [t-matrices/README.md](../../t-matrices/README.md)): first rows in `certificate.json`.
- Row sums of the four sequences: −53, 15, −117, 13 ((−53)² + 15² + (−117)² + 13² = 16892).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 163 103
t-matrices/T103.txt` writes the certificate.

## Status

**Known.** SageMath's development branch constructs this order since October 3, 2026 (merge of gh-42876, which added
the T-matrices of order 103 decoded from Schneider's matrix of order 2060; 11.0.beta1): `hadamard_matrix(16892)` uses
Cooper–Wallis with those T-matrices and the Williamson matrices of order 41. The matrix here is independent of it: it
comes from Paley's matrix of order 164 and our T-matrices of order 103, which were found separately and differ from
SageMath's.

In Table 1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 4223 has t = 3: the smallest known
Hadamard matrix with odd part 4223 had order 2³ · 4223. The order lies above the range of Cati–Pasechnik Table 4,
which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- Paley's constructions do not apply: 4n − 1 = 16891 = 7 · 19 · 127 and 2n − 1 = 8445 = 3 · 5 · 563 are not prime
  powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4); Yamada–Kiyasu does not apply: n − 2 = 4221 = 3² · 7 · 67 is not a prime power.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 4223 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 16892
```

SHA-256 of `H16892.csv`: `96a68f6a418f9bb41a736837a7ccc8c6af40898cd5a694d858394b8d0ee44b8f`

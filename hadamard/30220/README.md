# Hadamard matrix of order 30220

30220 = 4 · 7555, with 7555 = 5 · 1511.

## Construction

Paley's skew-Hadamard matrix of order 6044 (q = 6043, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 1511 once its points are listed along the four orbits of a torus element of PSL(2, 6043) and rescaled by
signs. Multiplying it by T-matrices of order 5 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 7555 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: the square of x ↦ −1/(x + τ) with τ = 3; sign normalisation c = +1.
- T-matrices of order 5 (from base sequences BS(3, 2), built by the program): first rows `++000`, `00+00`, `000+-`,
  `00000`.
- Row sums of the four sequences: 53, −101, −131, −7 (53² + (−101)² + (−131)² + (−7)² = 30220).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 6043 5` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 7555 has t = 3: the smallest known Hadamard matrix with
odd part 7555 had order 2³ · 7555. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 6044 is equivalent
  to an array of sixteen circulant blocks of order 1511 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and
  Yamada (1990, Theorem 4) multiply such an array by T-matrices of order 5.
- Paley's constructions do not apply: 4n − 1 = 30219 = 3 · 7 · 1439 and 2n − 1 = 15109 = 29 · 521 are not prime
  powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4); Yamada–Kiyasu does not apply: n − 2 = 7553 = 7 · 13 · 83 is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(30220, existence=True)`
  returns Unknown.

We found no construction or statement of order 30220. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 7555 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 30220
```

SHA-256 of `H30220.csv`: `1b03d02b8076dd7ce58ef7d528e948e797e5a2b923c0bf78055dc631dbf51173`

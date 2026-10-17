# Hadamard matrix of order 20260

20260 = 4 · 5065, with 5065 = 5 · 1013.

## Construction

Paley's skew-Hadamard matrix of order 4052 (q = 4051, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 1013 once its points are listed along the four orbits of a torus element of PSL(2, 4051) and rescaled by
signs. Multiplying it by T-matrices of order 5 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 5065 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = +1.
- T-matrices of order 5 (from base sequences BS(3, 2), built by the program): first rows `++000`, `00+00`, `000+-`,
  `00000`.
- Row sums of the four sequences: −43, 91, 43, −91 ((−43)² + 91² + 43² + (−91)² = 20260).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 4051 5` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 5065 has t = 3: the smallest known Hadamard matrix with
odd part 5065 had order 2³ · 5065. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 4052 is equivalent
  to an array of sixteen circulant blocks of order 1013 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and
  Yamada (1990, Theorem 4) multiply such an array by T-matrices of order 5.
- Paley's constructions do not apply: 4n − 1 = 20259 = 3² · 2251 and 2n − 1 = 10129 = 7 · 1447 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1 = 5066, and none is known to
  us (5065 is not a prime power); Yamada–Kiyasu does not apply: n − 2 = 5063 = 61 · 83 is not a prime power.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(20260, existence=True)`
  raises UnknownError in its Miyamoto test (no symmetric conference matrix of order 5066 is known to it), and returns
  Unknown when that Unknown is treated as False.

We found no construction or statement of order 20260. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 5065 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 20260
```

SHA-256 of `H20260.csv`: `32c154ea25e44094d7b0e77000184f9146a5bf68acb09c291692733e5b08c05d`

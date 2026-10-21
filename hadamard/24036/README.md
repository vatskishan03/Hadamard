# Hadamard matrix of order 24036

24036 = 4 · 6009, with 6009 = 3 · 2003.

## Construction

Paley's skew-Hadamard matrix of order 8012 (q = 8011, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 2003 once its points are listed along the four orbits of a torus element of PSL(2, 8011) and rescaled by
signs. Multiplying it by T-matrices of order 3 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 6009 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = +1.
- T-matrices of order 3 (from base sequences BS(2, 1), built by the program): first rows `+00`, `0+0`, `00+`, `000`.
- Row sums of the four sequences: 105, −103, 1, 49 (105² + (−103)² + 1² + 49² = 24036).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 8011 3` writes the
certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 6009 has t = 7: the smallest known Hadamard matrix with
odd part 6009 had order 2⁷ · 6009. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 8012 is equivalent
  to an array of sixteen circulant blocks of order 2003 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and
  Yamada (1990, Theorem 4) multiply such an array by T-matrices of order 3.
- Paley's constructions do not apply: 4n − 1 = 24035 = 5 · 11 · 19 · 23 and 2n − 1 = 12017 = 61 · 197 are not prime
  powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 6010 cannot exist (6009 is not a sum
  of two squares); Yamada–Kiyasu does not apply: n − 2 = 6007 is ≡ 7 (mod 8).
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (17 October 2026) does not construct it: `hadamard_matrix(24036, existence=True)`
  returns Unknown.

We found no construction or statement of order 24036. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 6009 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 24036
```

SHA-256 of `H24036.csv`: `d9b452032deaf2838936d5a56ff0ae0a58fa81536a70d2c6e06f430cb0e634f4`

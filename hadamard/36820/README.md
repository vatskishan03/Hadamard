# Hadamard matrix of order 36820

36820 = 4 · 9205, with 9205 = 5 · 7 · 263.

## Construction

Paley's skew-Hadamard matrix of order 1052 (q = 1051, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 263 once its points are listed along the four orbits of a torus element of PSL(2, 1051) and rescaled by
signs. Multiplying it by T-matrices of order 35 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 9205 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = −1.
- T-matrices of order 35 ([`t-matrices/T35.txt`](../../t-matrices/T35.txt), T-sequences of length 35 as stored in
  SageMath's `t_sequences` module): first rows `++++-++-++--00000000000000000000000`,
  `000000000000++++--+-+-+000000000000`, `00000000000000000000000+++0-0+0-00+`,
  `00000000000000000000000000+0-0-0-+0`.
- Row sums of the four sequences: 1, 35, 37, 185 (1² + 35² + 37² + 185² = 36820).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 1051 35
t-matrices/T35.txt` writes the certificate.

## Status

**New, beyond the tables, but obtainable classically** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#status)). In Table 1.53
of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 9205 has t = 3: the smallest known Hadamard matrix with
odd part 9205 had order 2³ · 9205. The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999.

- The entry does not appear in Đoković's updates (arXiv:0912.5091, arXiv:1008.2043).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 1052 is equivalent
  to an array of sixteen circulant blocks of order 263 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and Yamada
  (1990, Theorem 4) multiply such an array by T-matrices of order 35.
- Paley's constructions do not apply: 4n − 1 = 36819 = 3² · 4091 and 2n − 1 = 18409 = 41 · 449 are not prime powers.
- Miyamoto's theorem does not apply: a symmetric conference matrix of order n + 1 = 9206 cannot exist (9205 is not a sum
  of two squares); Yamada–Kiyasu does not apply: n − 2 = 9203 is ≡ 3 (mod 8).
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(36820, existence=True)`
  returns Unknown.

We found no construction or statement of order 36820. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 9205 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 36820
```

SHA-256 of `H36820.csv`: `5906a69099dce65f95dd3a6dbb907bc4d7c8857018a2a6ae4891eb1a4f5d0857`

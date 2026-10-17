# Hadamard matrix of order 31316

31316 = 4 · 7829, with 7829 prime and 7829 ≡ 1 (mod 4).

## Construction

Miyamoto's theorem with q = 7829, as described for [order 21748](../21748). The inputs are

- Paley's symmetric conference matrix of order 7830 over ℤ₇₈₂₉;
- a Hadamard matrix K of order 7828 = 4 · 103 · 19. It is the Cooper–Wallis product of the T-matrices of order 103 in
  [`t-matrices/T103.txt`](../../t-matrices/T103.txt), the ones used for [order 2060](../2060), and the Williamson
  matrices of order 19 listed in SageMath.

`certificate.json` contains q, the choice of conference matrix, the complete certificate of K and the SHA-256 of the
result.

## Status

**Known.** SageMath's development branch constructs this order since October 3, 2026 (merge of gh-42876, which added
the T-matrices of order 103 decoded from Schneider's matrix of order 2060; 11.0.beta1): `hadamard_matrix(31316)`
reaches Miyamoto's theorem with q = 7829 and a Cooper–Wallis matrix of order 7828 from those T-matrices and the
Williamson matrices of order 19. That predates our commit (October 12); we missed it, and the novelty re-audit of
October 16 found it. The matrix here is independent of it: our T-matrices of order 103 were found separately and differ
from SageMath's.

The Handbook (2007, Table 1.53) gives n = 7829 as t = 8: the smallest known Hadamard
matrix with odd part 7829 had order 2⁸ · 7829. The order lies above Cati–Pasechnik Table 4, which stops at n = 2999.

- The input order 4 · 1957 has t = 3 in the Handbook and in Cati–Pasechnik Table 4.
- No classical family gives order 4 · 7829.
- The order has followed from Miyamoto's theorem since JD-Jones-ASES/Hadamard-T published a matrix of order 7828 on
  August 29, 2026, but we found no construction or statement of it.
- T-matrices of order 103 first became public in 2026: JD-Jones-ASES/Hadamard-T decoded them on August 29, 2026, from
  the Hadamard matrix of order 2060 that Schneider posted. Ours were found independently (see [order 2060](../2060)).

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, conference matrix, certificate of K, SHA-256 |

The matrix is not stored. Its CSV encoding takes 2.5 GB (152 MB gzipped).

## Verification

```
python3 verify/verify.py --write 31316
```

This rebuilds and checks the matrix and writes `hadamard/31316/H31316.csv`. It takes 80 seconds on six cores and 2 GB
of memory.

SHA-256 of `H31316.csv`: `2ec73b0a165d4024ff33420418f2c413d5d2e7fceedf1d5e793f2d2dc7102c10`

# Hadamard matrix of order 32316

32316 = 4 · 8079, with 8079 = 3 · 2693.

## Construction

Paley's skew-Hadamard matrix of order 10772 (q = 10771, a prime ≡ 3 mod 8) is a 4 × 4 array of circulant blocks of
order s = 2693 once its points are listed along the four orbits of a torus element of PSL(2, 10771) and rescaled by
signs. Multiplying it by T-matrices of order 3 (Seberry–Yamada 1990, Theorem 4) gives four circulant ±1 sequences of
length 8079 with Σ PAF = 0 off 0, and the Goethals–Seidel array gives the Hadamard matrix. The construction and its
proof are in [PALEY_BLOCKS.md](../PALEY_BLOCKS.md).

- Torus element: x ↦ −1/(x + τ) with τ = 4; sign normalisation c = −1.
- T-matrices of order 3 (from base sequences BS(2, 1), built by the program): first rows `+00`, `0+0`, `00+`, `000`.
- Row sums of the four sequences: 113, −93, −17, 103 (113² + (−93)² + (−17)² + 103² = 32316).

## Search

None: the matrix is determined by q and the T-matrices. `python3 search/paley_blocks.py certify 10771 3` writes the
certificate.

## Status

**First explicit construction** ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md#orders-recorded-only-through-corollary-88i)).
In Table 1.53 of the Handbook of Combinatorial Designs (2nd ed., 2007), n = 8079 has t = 3, and Đoković's update
arXiv:1008.2043 lists it as settled, through Seberry–Yamada's Corollary 8.8(i) only (below). The order lies above the
range of Cati–Pasechnik Table 4, which stops at n = 2999. Its ingredients were all known before 2026.

- Đoković's update arXiv:1008.2043 lists it (Proposition 3.3), but only through Seberry–Yamada's Corollary 8.8(i),
  whose printed proof does not give Williamson-type matrices (see below).
- It follows from two published theorems that we found combined nowhere: Paley's matrix of order 10772 is equivalent
  to an array of sixteen circulant blocks of order 2693 (Schmidt 1999, after Ito; Yamada 1991), and Seberry and
  Yamada (1990, Theorem 4) multiply such an array by T-matrices of order 3.
- Paley's constructions do not apply: 4n − 1 = 32315 = 5 · 23 · 281 and 2n − 1 = 16157 = 107 · 151 are not prime
  powers.
- Miyamoto and Yamada–Kiyasu do not apply: 8079 is not a prime power, and neither is n − 2 = 8077 = 41 · 197.
- Cooper–Wallis does not apply: for each factorization n = a · b it needs T-matrices of order a and Williamson-type
  matrices of order b, and no such pair is known to us.
- SageMath's development version (October 17, 2026) does not construct it: `hadamard_matrix(32316, existence=True)`
  returns Unknown.

Đoković (arXiv:1008.2043, Proposition 3.3) obtained 8079 = 3 · 2693 from Williamson-type matrices of order 2693, which
he took from Corollary 8.8(i) of Seberry and Yamada (1992; Corollary 29(i) of JCMCC 7, 1990) with Williamson-type
matrices of order 673 = (2693 − 1)/4. The printed proof of that corollary does not give pairwise amicable matrices for
any Williamson-type input (in particular for Williamson-type matrices of order 673), and we found no other proof or
construction of Williamson-type matrices of order 2693 ([NOTES.md](../../NOTES.md) §7). Miyamoto's 1991 paper, credited
with this case, was not accessible to us. We found no valid construction of order 32316 before this one. The sources
searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power, c, the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences of length 8079 |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 32316
```

SHA-256 of `H32316.csv`: `a80be5280c6771adce3594d7b2fe2ea80d5bf444445e13bb97eb9eceb0adb918`

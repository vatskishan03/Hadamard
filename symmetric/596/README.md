# Symmetric Hadamard matrix of order 596

596 = 4 · 149, with 149 a prime ≡ 5 (mod 8).

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is the case q = 149 of [PERIODIC_GOLAY.md](../PERIODIC_GOLAY.md).
The construction takes Paley's conference matrix of order 150 and a Hadamard matrix K = [[P, Q], [Q, P]] of order 148
from a periodic Golay pair of length 74. H has an involutory automorphism with 296 fixed points.

- **The input.** Đoković and Kotsireas's periodic Golay pair of length 74 (Numer. Algorithms 69 (2015) 523–530,
  arXiv:1310.5773, parameter set (74; 36, 31; 30)). The two sets are unions of orbits of H = {1, 47, 63} on ℤ₇₄, with
  representatives {1, 4, 6, 7, 9, 12, 22, 23, 28, 29, 34, 42} and {1, 2, 4, 6, 9, 12, 17, 21, 22, 37, 55}. K is
  [[A, B], [−Bᵀ, Aᵀ]] with A, B circulant, reordered by the cosets of ⟨37⟩.
- **The field.** GF(149) = ℤ₁₄₉, χ the quadratic character, ν = 2 (the least non-square).
- **The Q-part.** On the 148 points x ≠ 0: A[x, y] = χ(x + y), Z[x, y] = [x + y = 0],
  F = [[A + Z, A − Z], [A − Z, A + Z]], B = [(𝟙; 𝟙) | (χ; χ) | (k_j; −k_j)] with k_j the columns of K (rows
  relabelled: row r is the point r + 1, row 74 + r the point −(r + 1)), and P = −BᵀFB/592.
- **The V-part.** V₀[x, y] = χ(xy − 2), V₀[x, ∞] = χ(x), V₀[∞, ∞] = 1, relabelled so that its zeros are the support
  of P.
- **The matrix.** H = [[F, B, B], [Bᵀ, P + V, P − V], [Bᵀ, P − V, P + V]].

## Search

None: the matrix is determined by q and the published pair. The pair is stored in `certificate.json` and hard-coded in
`search/sym_involution.py` (`PAIRS[149]`), and `python3 search/sym_involution.py certify 149` builds the certificate
from it.

## Status

**New** as a symmetric Hadamard matrix.

- The Đoković–Seberry census (2025-05-28) lists 149 among its undecided v for symmetric matrices.
  [PROPUS_WEIGHING.md](../PROPUS_WEIGHING.md) needs q ≡ 1 (mod 8) and does not apply to q = 149.
- The input is published (Đoković and Kotsireas); the construction from a non-symmetric K of the form [[P, Q], [Q, P]]
  in [PERIODIC_GOLAY.md](../PERIODIC_GOLAY.md#status) is new.
- Ordinary Hadamard matrices of order 596 are classical. The symmetry is what is new.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, the group, t, the periodic Golay pair (a, b), and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/596
```

`verify/verify.py` (construction `symmetric-involution`) checks the pair and K, builds the Q-part and V and checks
their identities, assembles H and checks H = Hᵀ, H Hᵀ = 596 I and the SHA-256 of the canonical CSV.

SHA-256 of `H596.csv`: `e7297eb866b7295bfee0bbfef3fd8e3b292fec67e7cbdadc975d5003db5b35c5`

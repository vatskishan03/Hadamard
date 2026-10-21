# Symmetric Hadamard matrix of order 404

404 = 4 · 101, with 101 a prime ≡ 5 (mod 8).

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is the case q = 101 of [PERIODIC_GOLAY.md](../PERIODIC_GOLAY.md).
The construction takes Paley's conference matrix of order 102 and a Hadamard matrix K = [[P, Q], [Q, P]] of order 100
from a periodic Golay pair of length 50. H has an involutory automorphism with 200 fixed points.

- **The input.** A periodic Golay pair (a, b) of length 50 over Γ = ℤ₂ × ℤ₅² (PAF_a + PAF_b = 100 δ), given in
  `certificate.json` in the mixed-radix order of Γ. It comes from our search over GF(25) with A, B, C, D constant on
  the classes of index 8, and its K has t = (1, 0, 0). K is [[A, B], [−Bᵀ, Aᵀ]] with A, B group-developed over Γ,
  reordered by the cosets of ⟨t⟩.
- **The field.** GF(101) = ℤ₁₀₁, χ the quadratic character, ν = 2 (the least non-square).
- **The Q-part.** On the 100 points x ≠ 0: A[x, y] = χ(x + y), Z[x, y] = [x + y = 0],
  F = [[A + Z, A − Z], [A − Z, A + Z]], B = [(𝟙; 𝟙) | (χ; χ) | (k_j; −k_j)] with k_j the columns of K (rows
  relabelled: row r is the point r + 1, row 50 + r the point −(r + 1)), and P = −BᵀFB/400.
- **The V-part.** V₀[x, y] = χ(xy − 2), V₀[x, ∞] = χ(x), V₀[∞, ∞] = 1, relabelled so that its zeros are the support
  of P.
- **The matrix.** H = [[F, B, B], [Bᵀ, P + V, P − V], [Bᵀ, P − V, P + V]].

## Search

The periodic Golay pair of length 50 over ℤ₂ × ℤ₅² comes from an exhaustive enumeration of ternary quadruples A, B, C,
D constant on the cosets of the index-8 subgroup of GF(25)* with complementary supports (432128 pairs after
power-spectrum pruning). The first pair was converted to (a, b), and the matrix is then determined.

The pair is stored in `certificate.json` and hard-coded in `search/sym_involution.py` (`PAIRS[101]`), and
`python3 search/sym_involution.py certify 101` builds the certificate from it. The certificate alone proves the
result; the search programs for this order are not included in this repository.

## Status

**New** as a symmetric Hadamard matrix.

- The Đoković–Seberry census (2025-05-28) lists 101 among its undecided v for symmetric matrices, and it is one of its
  smallest open values. [PROPUS_WEIGHING.md](../PROPUS_WEIGHING.md) needs q ≡ 1 (mod 8) and does not apply to q = 101.
- The construction from a non-symmetric K of the form [[P, Q], [Q, P]] is new
  ([PERIODIC_GOLAY.md](../PERIODIC_GOLAY.md#status)). Periodic Golay pairs of length 50 were found earlier by Đoković
  (Ann. Comb. 15 (2011) 233–254) and by Kotsireas and Koukouvinos (Int. J. Appl. Math. 21 (2008) 509–514), as listed
  in arXiv:1310.5773; the pair used here comes from our search. A published cyclic pair of length 50 would give K as
  well.
- Ordinary Hadamard matrices of order 404 are classical. The symmetry is what is new.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, the group, t, the periodic Golay pair (a, b), and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/404
```

`verify/verify.py` (construction `symmetric-involution`) checks the pair and K, builds the Q-part and V and checks
their identities, assembles H and checks H = Hᵀ, H Hᵀ = 404 I and the SHA-256 of the canonical CSV.

SHA-256 of `H404.csv`: `291df30dc0701aed01e3ab82e9ee78c1145a3bb7e8286595c151eb39d4422dcf`

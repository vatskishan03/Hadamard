# Symmetric Hadamard matrix of order 6684

6684 = 4 · 1671, with 1671 = 3 · 557; 1671 = q + 2 with q = 1669 a prime ≡ 5 (mod 8) and (q + 3)/2 = 836.

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is the case q = 1669 of the construction in
[YAMADA_KIYASU.md](../YAMADA_KIYASU.md): Yamada's (Kiyasu's) Lemma 26 data put in quaternary form, which is symmetric for
any skew input, and Turyn's map to a real matrix.

- **The skew input.** Our skew-Hadamard matrix of order 836 ([skew/836](../../skew/836)). Normalised to first row all +1, its core S of order 835 has
  Sᵀ = −S, Se = 0 and S² = J − 835I.
- **Yamada's circulants.** M and N of order 835, with first rows β and α of Yamada's Lemma 3 over
  GF(1669²) = GF(1669)[w]/(w² − 2) and ξ = 2 + w primitive: symmetric, M of zero diagonal,
  M² + N² = 1669I.
- **The quaternary matrix.** With G = ½((M + N) + i(M − N)), P = [[1, −1], [−1, 1]] and e the all-one row of length 1670,
  L = [[−1, −e, −1, ie], [−eᵀ, ½(1 + i)I ⊗ J₂ + G ⊗ P, −ieᵀ, −iI ⊗ P − S ⊗ J₂], [−1, −ie, 1, e],
  [ieᵀ, −iI ⊗ P + S ⊗ J₂, eᵀ, ½(i − 1)I ⊗ J₂ + Ḡ ⊗ P]] of order 3342 is symmetric with entries ±1, ±i and
  L L̄ = 3342I.
- **The matrix.** With L = X + iY, H = [[X + Y, X − Y], [X − Y, −X − Y]].

## Search

None: the matrix is determined by q, r, ξ and the skew input. `python3 search/sym_yk.py certify 1671` writes the
certificate.

## Status

**New** as a symmetric Hadamard matrix.

- The Đoković–Seberry census (2025-05-28) lists 1671 among its 595 undecided v for symmetric matrices. Seberry–Yamada's
  Corollary 27(i) in the propus array would need a symmetric conference matrix of order 836 or a skew-Hadamard
  matrix of order 836 with an involutive anti-automorphism ([NOTES.md](../../NOTES.md) §29).
- The symmetric quaternary form of Yamada's matrix is the new step ([YAMADA_KIYASU.md](../YAMADA_KIYASU.md#status)); we
  found it in none of the sources listed there. The skew input is ours (2026).
- Ordinary Hadamard matrices of order 6684 are classical (Yamada 1989, Theorem 4, from the same skew input); the
  symmetry is what is new.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, r, ξ, the complete certificate of the skew input, and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/6684
```

`verify/verify.py` (construction `symmetric-yamada-kiyasu`) rebuilds the skew input from its certificate and checks it
(Q + Qᵀ = 2I), computes Yamada's circulants and checks M² + N² = 1669I, builds X and Y and checks that they are symmetric
with disjoint supports covering every position, and checks H = Hᵀ, H Hᵀ = 6684 I and the SHA-256 of the canonical CSV.

SHA-256 of `H6684.csv`: `86688903aae69519fb81b8066aca5e687037f7c8b9a14e29b056e79fcd46e8fe`

# Symmetric Hadamard matrix of order 3772

3772 = 4 · 943, with 943 = 23 · 41; 943 = q + 2 with q = 941 a prime ≡ 5 (mod 8) and (q + 3)/2 = 472.

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is the case q = 941 of the construction in
[YAMADA_KIYASU.md](../YAMADA_KIYASU.md): Yamada's (Kiyasu's) Lemma 26 data put in quaternary form, which is symmetric for
any skew input, and Turyn's map to a real matrix.

- **The skew input.** Fletcher, Koukouvinos and Seberry's skew Goethals–Seidel matrix of order 236 (2004; as stored in SageMath), doubled once. Normalised to first row all +1, its core S of order 471 has
  Sᵀ = −S, Se = 0 and S² = J − 471I.
- **Yamada's circulants.** M and N of order 471, with first rows β and α of Yamada's Lemma 3 over
  GF(941²) = GF(941)[w]/(w² − 2) and ξ = 2 + w primitive: symmetric, M of zero diagonal,
  M² + N² = 941I.
- **The quaternary matrix.** With G = ½((M + N) + i(M − N)), P = [[1, −1], [−1, 1]] and e the all-one row of length 942,
  L = [[−1, −e, −1, ie], [−eᵀ, ½(1 + i)I ⊗ J₂ + G ⊗ P, −ieᵀ, −iI ⊗ P − S ⊗ J₂], [−1, −ie, 1, e],
  [ieᵀ, −iI ⊗ P + S ⊗ J₂, eᵀ, ½(i − 1)I ⊗ J₂ + Ḡ ⊗ P]] of order 1886 is symmetric with entries ±1, ±i and
  L L̄ = 1886I.
- **The matrix.** With L = X + iY, H = [[X + Y, X − Y], [X − Y, −X − Y]].

## Search

None: the matrix is determined by q, r, ξ and the skew input. `python3 search/sym_yk.py certify 943` writes the
certificate.

## Status

**Concurrent** as a symmetric Hadamard matrix; we claim no priority for the order. Zaitsev, Kharaghani,
Mohammadian and Tayfeh-Rezaie's page of quaternary Hadamard matrices (last updated 14 October 2026) lists
QH1886 under "unpublished approaches". It is a quaternary Yamada–Kiyasu matrix of order 1886, a monomial
matrix takes it to a symmetric one (found with nauty, checked exactly: same rows, columns up to units, symmetric,
Hadamard), and Turyn's map then gives a symmetric Hadamard matrix of order 3772. The page does not mention
symmetry or real matrices, and the construction for an arbitrary skew input is ours
([YAMADA_KIYASU.md](../YAMADA_KIYASU.md#status)).

- The Đoković–Seberry census (2025-05-28) lists 943 among its 595 undecided v for symmetric matrices. Seberry–Yamada's
  Corollary 27(i) in the propus array would need a symmetric conference matrix of order 472 or a skew-Hadamard
  matrix of order 472 with an involutive anti-automorphism ([NOTES.md](../../NOTES.md) §29).
- The symmetric quaternary form of Yamada's matrix is the new step ([YAMADA_KIYASU.md](../YAMADA_KIYASU.md#status)); we
  found it in none of the sources listed there. The skew input is classical.
- Ordinary Hadamard matrices of order 3772 are classical (Yamada 1989, Theorem 4, from the same skew input); the
  symmetry is what is new.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, r, ξ, the complete certificate of the skew input, and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/3772
```

`verify/verify.py` (construction `symmetric-yamada-kiyasu`) rebuilds the skew input from its certificate and checks it
(Q + Qᵀ = 2I), computes Yamada's circulants and checks M² + N² = 941I, builds X and Y and checks that they are symmetric
with disjoint supports covering every position, and checks H = Hᵀ, H Hᵀ = 3772 I and the SHA-256 of the canonical CSV.

SHA-256 of `H3772.csv`: `2a8bb1650133f259df9b73a7d20d906eaa4fd45a8c4731b15d531b289758c6af`

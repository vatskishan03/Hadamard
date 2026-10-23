# Symmetric Hadamard matrix of order 7612

7612 = 4 · 1903, with 1903 = 11 · 173; 1903 = q + 2 with q = 1901 a prime ≡ 5 (mod 8) and (q + 3)/2 = 952.

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is the case q = 1901 of the construction in
[YAMADA_KIYASU.md](../YAMADA_KIYASU.md): Yamada's (Kiyasu's) Lemma 26 data put in quaternary form, which is symmetric for
any skew input, and Turyn's map to a real matrix.

- **The skew input.** Our skew-Hadamard matrix of order 476 ([skew/476](../../skew/476)), doubled once
  ([[S, S], [−Sᵀ, Sᵀ]]) to order 952. Normalised to first row all +1, its core S of order 951 has Sᵀ = −S, Se = 0 and
  S² = J − 951I.
- **Yamada's circulants.** M and N of order 951, with first rows β and α of Yamada's Lemma 3 over
  GF(1901²) = GF(1901)[w]/(w² − 2) and ξ = 8 + w primitive: symmetric, M of zero diagonal, M² + N² = 1901I.
- **The quaternary matrix.** With G = ½((M + N) + i(M − N)), P = [[1, −1], [−1, 1]] and e the all-one row of length 1902,
  L = [[−1, −e, −1, ie], [−eᵀ, ½(1 + i)I ⊗ J₂ + G ⊗ P, −ieᵀ, −iI ⊗ P − S ⊗ J₂], [−1, −ie, 1, e],
  [ieᵀ, −iI ⊗ P + S ⊗ J₂, eᵀ, ½(i − 1)I ⊗ J₂ + Ḡ ⊗ P]] of order 3806 is symmetric with entries ±1, ±i and
  L L̄ = 3806I.
- **The matrix.** With L = X + iY, H = [[X + Y, X − Y], [X − Y, −X − Y]].

## Search

None: the matrix is determined by q, r, ξ and the skew input. `python3 search/sym_yk.py certify 1903` writes the
certificate.

## Status

**New** as a symmetric Hadamard matrix.

- The Đoković–Seberry census (2025-05-28) lists 1903 among its 595 undecided v for symmetric matrices. Seberry–Yamada's
  Corollary 27(i) in the propus array would need (q + 1)/2 = 951 = 3 · 317 to be a prime power, or a skew-Hadamard
  matrix of order 952 with an involutive anti-automorphism ([NOTES.md](../../NOTES.md) §29).
- The symmetric quaternary form of Yamada's matrix is the new step ([YAMADA_KIYASU.md](../YAMADA_KIYASU.md#status)); we
  found it in none of the sources listed there. The skew input is ours (2026): before it, no skew-Hadamard matrix of
  order 952 = 8 · 119 was known (Cati–Pasechnik Table 5 gives 119(4)).
- An ordinary Hadamard matrix of order 7612 is in [hadamard/7612](../../hadamard/7612) (Paley's matrix in circulant
  blocks times T-matrices, a first explicit construction), and Yamada's Theorem 4 now gives another from the same skew
  input; the symmetry is what is new.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, r, ξ, the complete certificate of the skew input, and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/7612
```

`verify/verify.py` (construction `symmetric-yamada-kiyasu`) rebuilds the skew input from its certificate and checks it
(Q + Qᵀ = 2I), computes Yamada's circulants and checks M² + N² = 1901I, builds X and Y and checks that they are symmetric
with disjoint supports covering every position, and checks H = Hᵀ, H Hᵀ = 7612 I and the SHA-256 of the canonical CSV.

SHA-256 of `H7612.csv`: `10df018b1bdc0b04867beee9509a5beb6ced8c33f3c6015c096c851e67cbe480`

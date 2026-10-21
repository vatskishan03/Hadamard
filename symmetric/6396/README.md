# Symmetric Hadamard matrix of order 6396

6396 = 4 · 1599, with 1599 = 3 · 13 · 41; 1599 = q + 2 with q = 1597 a prime ≡ 5 (mod 8) and (q + 3)/2 = 800.

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is the case q = 1597 of the construction in
[YAMADA_KIYASU.md](../YAMADA_KIYASU.md): Yamada's (Kiyasu's) Lemma 26 data put in quaternary form, which is symmetric
for any skew input, and Turyn's map to a real matrix.

- **The skew input.** Paley's skew-Hadamard matrix I + C of order 200 (p = 199), doubled twice ([[S, S], [−Sᵀ, Sᵀ]]).
  Normalised to first row all +1, its core S of order 799 has Sᵀ = −S, Se = 0 and S² = J − 799I.
- **Yamada's circulants.** M and N of order 799, with first rows β and α of Yamada's Lemma 3 over
  GF(1597²) = GF(1597)[w]/(w² − 2) and ξ = 6 + w primitive: symmetric, M of zero diagonal,
  M² + N² = 1597I.
- **The quaternary matrix.** With G = ½((M + N) + i(M − N)), P = [[1, −1], [−1, 1]] and e the all-one row of length
  1598, L = [[−1, −e, −1, ie], [−eᵀ, ½(1 + i)I ⊗ J₂ + G ⊗ P, −ieᵀ, −iI ⊗ P − S ⊗ J₂], [−1, −ie, 1, e],
  [ieᵀ, −iI ⊗ P + S ⊗ J₂, eᵀ, ½(i − 1)I ⊗ J₂ + Ḡ ⊗ P]] of order 3198 is symmetric with entries ±1, ±i and
  L L̄ = 3198I.
- **The matrix.** With L = X + iY, H = [[X + Y, X − Y], [X − Y, −X − Y]].

## Search

None: the matrix is determined by q, r, ξ and the skew input.

**Reproduction.**

| Command | Output |
|---|---|
| `python3 search/sym_yk.py certify 1599` | `symmetric/6396/certificate.json`, identical to the stored file (the program refuses to overwrite it: run it in a copy of the repository without that file) |
| `python3 verify/verify.py symmetric/6396` | `ok`, SHA-256 `95d940b9ce2d4404954de087a81d050567833486030783b07b1a87afd00a6454` |

## Status

**New** as a symmetric Hadamard matrix.

- The Đoković–Seberry census (2025-05-28) lists 1599 among its 595 undecided v for symmetric matrices.
  Seberry–Yamada's Corollary 27(i) in the propus array would need a skew-Hadamard matrix of order 800 with an
  involutive anti-automorphism ([NOTES.md](../../NOTES.md) §29).
- The symmetric quaternary form of Yamada's matrix is the new step ([YAMADA_KIYASU.md](../YAMADA_KIYASU.md#status));
  we found it in none of the sources listed there. The skew input is classical.
- Ordinary Hadamard matrices of order 6396 are classical (Yamada 1989, Theorem 4, from the same skew input); the
  symmetry is what is new.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, r, ξ, the complete certificate of the skew input, and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/6396
```

`verify/verify.py` (construction `symmetric-yamada-kiyasu`) rebuilds the skew input from its certificate and checks it
(Q + Qᵀ = 2I), computes Yamada's circulants and checks M² + N² = 1597I, builds X and Y and checks that they are
symmetric with disjoint supports covering every position, and checks H = Hᵀ, H Hᵀ = 6396 I and the SHA-256 of the
canonical CSV.

SHA-256 of `H6396.csv`: `95d940b9ce2d4404954de087a81d050567833486030783b07b1a87afd00a6454`

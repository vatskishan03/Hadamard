# Symmetric Hadamard matrix of order 4924

4924 = 4 · 1231, with 1231 prime; 1231 = q + 2 with q = 1229 a prime ≡ 5 (mod 8) and (q + 3)/2 = 616.

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is the case q = 1229 of the construction in
[YAMADA_KIYASU.md](../YAMADA_KIYASU.md): Yamada's (Kiyasu's) Lemma 26 data put in quaternary form, which is symmetric
for any skew input, and Turyn's map to a real matrix.

- **The skew input.** Paley's skew-Hadamard matrix I + C of order 308 (p = 307), doubled once ([[S, S], [−Sᵀ, Sᵀ]]).
  Normalised to first row all +1, its core S of order 615 has Sᵀ = −S, Se = 0 and S² = J − 615I.
- **Yamada's circulants.** M and N of order 615, with first rows β and α of Yamada's Lemma 3 over
  GF(1229²) = GF(1229)[w]/(w² − 2) and ξ = 2 + w primitive: symmetric, M of zero diagonal,
  M² + N² = 1229I.
- **The quaternary matrix.** With G = ½((M + N) + i(M − N)), P = [[1, −1], [−1, 1]] and e the all-one row of length
  1230, L = [[−1, −e, −1, ie], [−eᵀ, ½(1 + i)I ⊗ J₂ + G ⊗ P, −ieᵀ, −iI ⊗ P − S ⊗ J₂], [−1, −ie, 1, e],
  [ieᵀ, −iI ⊗ P + S ⊗ J₂, eᵀ, ½(i − 1)I ⊗ J₂ + Ḡ ⊗ P]] of order 2462 is symmetric with entries ±1, ±i and
  L L̄ = 2462I.
- **The matrix.** With L = X + iY, H = [[X + Y, X − Y], [X − Y, −X − Y]].

## Search

None: the matrix is determined by q, r, ξ and the skew input.

**Reproduction.**

| Command | Output |
|---|---|
| `python3 search/sym_yk.py certify 1231` | `symmetric/4924/certificate.json`, identical to the stored file (the program refuses to overwrite it: run it in a copy of the repository without that file) |
| `python3 verify/verify.py symmetric/4924` | `ok`, SHA-256 `122d6d40599fd88226d1d19faca47faff95a089777b904ce6ba6ddd350509947` |

## Status

**New** as a symmetric Hadamard matrix.

- The Đoković–Seberry census (2025-05-28) lists 1231 among its 595 undecided v for symmetric matrices.
  Seberry–Yamada's Corollary 27(i) in the propus array would need a skew-Hadamard matrix of order 616 with an
  involutive anti-automorphism ([NOTES.md](../../NOTES.md) §29).
- The symmetric quaternary form of Yamada's matrix is the new step ([YAMADA_KIYASU.md](../YAMADA_KIYASU.md#status));
  we found it in none of the sources listed there. The skew input is classical.
- Ordinary Hadamard matrices of order 4924 are classical (Yamada 1989, Theorem 4, from the same skew input); the
  symmetry is what is new.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, r, ξ, the complete certificate of the skew input, and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/4924
```

`verify/verify.py` (construction `symmetric-yamada-kiyasu`) rebuilds the skew input from its certificate and checks it
(Q + Qᵀ = 2I), computes Yamada's circulants and checks M² + N² = 1229I, builds X and Y and checks that they are
symmetric with disjoint supports covering every position, and checks H = Hᵀ, H Hᵀ = 4924 I and the SHA-256 of the
canonical CSV.

SHA-256 of `H4924.csv`: `122d6d40599fd88226d1d19faca47faff95a089777b904ce6ba6ddd350509947`

# Symmetric Hadamard matrix of order 5980

5980 = 4 · 1495, with 1495 = 5 · 13 · 23; 1495 = q + 2 with q = 1493 a prime ≡ 5 (mod 8) and (q + 3)/2 = 748.

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is the case q = 1493 of the construction in
[YAMADA_KIYASU.md](../YAMADA_KIYASU.md): Yamada's (Kiyasu's) Lemma 26 data put in quaternary form, which is symmetric for
any skew input, and Turyn's map to a real matrix.

- **The skew input.** The two-block skew-Hadamard matrix of order 748 from Szekeres's (1969) complementary difference sets in GF(373) (fourth powers). Normalised to first row all +1, its core S of order 747 has
  Sᵀ = −S, Se = 0 and S² = J − 747I.
- **Yamada's circulants.** M and N of order 747, with first rows β and α of Yamada's Lemma 3 over
  GF(1493²) = GF(1493)[w]/(w² − 2) and ξ = 4 + w primitive: symmetric, M of zero diagonal,
  M² + N² = 1493I.
- **The quaternary matrix.** With G = ½((M + N) + i(M − N)), P = [[1, −1], [−1, 1]] and e the all-one row of length 1494,
  L = [[−1, −e, −1, ie], [−eᵀ, ½(1 + i)I ⊗ J₂ + G ⊗ P, −ieᵀ, −iI ⊗ P − S ⊗ J₂], [−1, −ie, 1, e],
  [ieᵀ, −iI ⊗ P + S ⊗ J₂, eᵀ, ½(i − 1)I ⊗ J₂ + Ḡ ⊗ P]] of order 2990 is symmetric with entries ±1, ±i and
  L L̄ = 2990I.
- **The matrix.** With L = X + iY, H = [[X + Y, X − Y], [X − Y, −X − Y]].

## Search

None: the matrix is determined by q, r, ξ and the skew input. `python3 search/sym_yk.py certify 1495` writes the
certificate.

## Status

**New** as a symmetric Hadamard matrix.

- The Đoković–Seberry census (2025-05-28) lists 1495 among its 595 undecided v for symmetric matrices. Seberry–Yamada's
  Corollary 27(i) in the propus array would need a symmetric conference matrix of order 748 or a skew-Hadamard
  matrix of order 748 with an involutive anti-automorphism ([NOTES.md](../../NOTES.md) §29).
- The symmetric quaternary form of Yamada's matrix is the new step ([YAMADA_KIYASU.md](../YAMADA_KIYASU.md#status)); we
  found it in none of the sources listed there. The skew input is classical.
- Ordinary Hadamard matrices of order 5980 are classical (Yamada 1989, Theorem 4, from the same skew input); the
  symmetry is what is new.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, r, ξ, the complete certificate of the skew input, and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/5980
```

`verify/verify.py` (construction `symmetric-yamada-kiyasu`) rebuilds the skew input from its certificate and checks it
(Q + Qᵀ = 2I), computes Yamada's circulants and checks M² + N² = 1493I, builds X and Y and checks that they are symmetric
with disjoint supports covering every position, and checks H = Hᵀ, H Hᵀ = 5980 I and the SHA-256 of the canonical CSV.

SHA-256 of `H5980.csv`: `b778f049c782fe666d26e10c64fd548d6b7b8720a9a7ce355b3c1650883f04c7`

# Symmetric Hadamard matrix of order 2716

2716 = 4 · 679, with 679 = 7 · 97; 679 = q + 2 with q = 677 a prime ≡ 5 (mod 8) and (q + 3)/2 = 340.

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is the case q = 677 of the construction in
[YAMADA_KIYASU.md](../YAMADA_KIYASU.md): Yamada's (Kiyasu's) Lemma 26 data put in quaternary form, which is symmetric
for any skew input, and Turyn's map to a real matrix.

- **The skew input.** The two-block skew-Hadamard matrix of order 340 from Whiteman's (1971) complementary difference
  sets in GF(13²) (eighth powers). Normalised to first row all +1, its core S of order 339 has Sᵀ = −S, Se = 0 and
  S² = J − 339I.
- **Yamada's circulants.** M and N of order 339, with first rows β and α of Yamada's Lemma 3 over
  GF(677²) = GF(677)[w]/(w² − 2) and ξ = 3 + w primitive: symmetric, M of zero diagonal,
  M² + N² = 677I.
- **The quaternary matrix.** With G = ½((M + N) + i(M − N)), P = [[1, −1], [−1, 1]] and e the all-one row of length
  678, L = [[−1, −e, −1, ie], [−eᵀ, ½(1 + i)I ⊗ J₂ + G ⊗ P, −ieᵀ, −iI ⊗ P − S ⊗ J₂], [−1, −ie, 1, e],
  [ieᵀ, −iI ⊗ P + S ⊗ J₂, eᵀ, ½(i − 1)I ⊗ J₂ + Ḡ ⊗ P]] of order 1358 is symmetric with entries ±1, ±i and
  L L̄ = 1358I.
- **The matrix.** With L = X + iY, H = [[X + Y, X − Y], [X − Y, −X − Y]].

## Search

None: the matrix is determined by q, r, ξ and the skew input.

**Reproduction.**

| Command | Output |
|---|---|
| `python3 search/sym_yk.py certify 679` | `symmetric/2716/certificate.json`, identical to the stored file (the program refuses to overwrite it: run it in a copy of the repository without that file) |
| `python3 verify/verify.py symmetric/2716` | `ok`, SHA-256 `878175fb3c0e5d5fa1a9e7dfbdbc5343c039cc2d1ba5dfb76645534276cc03f9` |

## Status

**Concurrent** as a symmetric Hadamard matrix; we claim no priority for the order. Zaitsev, Kharaghani, Mohammadian
and Tayfeh-Rezaie's page of quaternary Hadamard matrices (last updated 14 October 2026, unchanged on 22 October) lists
QH1358 under "unpublished approaches". It is a quaternary Yamada–Kiyasu matrix of order 1358, a monomial matrix takes
it to a symmetric one (found with nauty and checked exactly: same rows, columns up to units, symmetric, Hadamard; the
programs for this check are not included in this repository), and Turyn's map then gives a symmetric Hadamard matrix
of order 2716. The page does not mention symmetry or real matrices, and the construction for an arbitrary skew input
is ours ([YAMADA_KIYASU.md](../YAMADA_KIYASU.md#status)).

- The Đoković–Seberry census (2025-05-28) lists 679 among its 595 undecided v for symmetric matrices.
  Seberry–Yamada's Corollary 27(i) in the propus array would need a skew-Hadamard matrix of order 340 with an
  involutive anti-automorphism ([NOTES.md](../../NOTES.md) §29).
- The symmetric quaternary form of Yamada's matrix is the new step ([YAMADA_KIYASU.md](../YAMADA_KIYASU.md#status));
  we found it in none of the sources listed there. The skew input is classical.
- Ordinary Hadamard matrices of order 2716 are classical (Yamada 1989, Theorem 4, from the same skew input); the
  symmetry is what is new.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, r, ξ, the complete certificate of the skew input, and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/2716
```

`verify/verify.py` (construction `symmetric-yamada-kiyasu`) rebuilds the skew input from its certificate and checks it
(Q + Qᵀ = 2I), computes Yamada's circulants and checks M² + N² = 677I, builds X and Y and checks that they are
symmetric with disjoint supports covering every position, and checks H = Hᵀ, H Hᵀ = 2716 I and the SHA-256 of the
canonical CSV.

SHA-256 of `H2716.csv`: `878175fb3c0e5d5fa1a9e7dfbdbc5343c039cc2d1ba5dfb76645534276cc03f9`

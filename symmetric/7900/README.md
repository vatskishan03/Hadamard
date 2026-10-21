# Symmetric Hadamard matrix of order 7900

7900 = 4 · 1975, with 1975 = 5² · 79; 1975 = q + 2 with q = 1973 a prime ≡ 5 (mod 8) and (q + 3)/2 = 988.

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is the case q = 1973 of the construction in
[YAMADA_KIYASU.md](../YAMADA_KIYASU.md): Yamada's (Kiyasu's) Lemma 26 data put in quaternary form, which is symmetric
for any skew input, and Turyn's map to a real matrix.

- **The skew input.** Đoković's skew supplementary difference sets over ℤ₂₄₇ (2008; as stored in SageMath), a skew
  Goethals–Seidel matrix of order 988. Normalised to first row all +1, its core S of order 987 has Sᵀ = −S, Se = 0 and
  S² = J − 987I.
- **Yamada's circulants.** M and N of order 987, with first rows β and α of Yamada's Lemma 3 over
  GF(1973²) = GF(1973)[w]/(w² − 2) and ξ = 2 + w primitive: symmetric, M of zero diagonal,
  M² + N² = 1973I.
- **The quaternary matrix.** With G = ½((M + N) + i(M − N)), P = [[1, −1], [−1, 1]] and e the all-one row of length
  1974, L = [[−1, −e, −1, ie], [−eᵀ, ½(1 + i)I ⊗ J₂ + G ⊗ P, −ieᵀ, −iI ⊗ P − S ⊗ J₂], [−1, −ie, 1, e],
  [ieᵀ, −iI ⊗ P + S ⊗ J₂, eᵀ, ½(i − 1)I ⊗ J₂ + Ḡ ⊗ P]] of order 3950 is symmetric with entries ±1, ±i and
  L L̄ = 3950I.
- **The matrix.** With L = X + iY, H = [[X + Y, X − Y], [X − Y, −X − Y]].

## Search

None: the matrix is determined by q, r, ξ and the skew input.

**Reproduction.**

| Command | Output |
|---|---|
| `python3 search/sym_yk.py certify 1975` | `symmetric/7900/certificate.json`, identical to the stored file (the program refuses to overwrite it: run it in a copy of the repository without that file) |
| `python3 verify/verify.py symmetric/7900` | `ok`, SHA-256 `f8bb5233f1f0d39ca4024d9a285077fd26c2f5dde87ff03b4aa18973ca349f67` |

## Status

**New** as a symmetric Hadamard matrix.

- The Đoković–Seberry census (2025-05-28) lists 1975 among its 595 undecided v for symmetric matrices.
  Seberry–Yamada's Corollary 27(i) in the propus array would need a skew-Hadamard matrix of order 988 with an
  involutive anti-automorphism ([NOTES.md](../../NOTES.md) §29).
- The symmetric quaternary form of Yamada's matrix is the new step ([YAMADA_KIYASU.md](../YAMADA_KIYASU.md#status));
  we found it in none of the sources listed there. The skew input is classical.
- Ordinary Hadamard matrices of order 7900 are classical (Yamada 1989, Theorem 4, from the same skew input); the
  symmetry is what is new.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, r, ξ, the complete certificate of the skew input, and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/7900
```

`verify/verify.py` (construction `symmetric-yamada-kiyasu`) rebuilds the skew input from its certificate and checks it
(Q + Qᵀ = 2I), computes Yamada's circulants and checks M² + N² = 1973I, builds X and Y and checks that they are
symmetric with disjoint supports covering every position, and checks H = Hᵀ, H Hᵀ = 7900 I and the SHA-256 of the
canonical CSV.

SHA-256 of `H7900.csv`: `f8bb5233f1f0d39ca4024d9a285077fd26c2f5dde87ff03b4aa18973ca349f67`

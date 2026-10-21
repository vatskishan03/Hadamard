# Symmetric Hadamard matrix of order 7804

7804 = 4 · 1951, with 1951 prime; 1951 = q + 2 with q = 1949 a prime ≡ 5 (mod 8) and (q + 3)/2 = 976.

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is the case q = 1949 of the construction in
[YAMADA_KIYASU.md](../YAMADA_KIYASU.md): Yamada's (Kiyasu's) Lemma 26 data put in quaternary form, which is symmetric
for any skew input, and Turyn's map to a real matrix.

- **The skew input.** Paley's skew-Hadamard matrix I + C of order 488 (p = 487), doubled once ([[S, S], [−Sᵀ, Sᵀ]]).
  Normalised to first row all +1, its core S of order 975 has Sᵀ = −S, Se = 0 and S² = J − 975I.
- **Yamada's circulants.** M and N of order 975, with first rows β and α of Yamada's Lemma 3 over
  GF(1949²) = GF(1949)[w]/(w² − 2) and ξ = 3 + w primitive: symmetric, M of zero diagonal,
  M² + N² = 1949I.
- **The quaternary matrix.** With G = ½((M + N) + i(M − N)), P = [[1, −1], [−1, 1]] and e the all-one row of length
  1950, L = [[−1, −e, −1, ie], [−eᵀ, ½(1 + i)I ⊗ J₂ + G ⊗ P, −ieᵀ, −iI ⊗ P − S ⊗ J₂], [−1, −ie, 1, e],
  [ieᵀ, −iI ⊗ P + S ⊗ J₂, eᵀ, ½(i − 1)I ⊗ J₂ + Ḡ ⊗ P]] of order 3902 is symmetric with entries ±1, ±i and
  L L̄ = 3902I.
- **The matrix.** With L = X + iY, H = [[X + Y, X − Y], [X − Y, −X − Y]].

## Search

None: the matrix is determined by q, r, ξ and the skew input.

**Reproduction.**

| Command | Output |
|---|---|
| `python3 search/sym_yk.py certify 1951` | `symmetric/7804/certificate.json`, identical to the stored file (the program refuses to overwrite it: run it in a copy of the repository without that file) |
| `python3 verify/verify.py symmetric/7804` | `ok`, SHA-256 `c840b35bbe34860e218542a9ad54136d3c6aa16d5078af0e978e03f8f3062096` |

## Status

**New** as a symmetric Hadamard matrix.

- The Đoković–Seberry census (2025-05-28) lists 1951 among its 595 undecided v for symmetric matrices.
  Seberry–Yamada's Corollary 27(i) in the propus array would need a skew-Hadamard matrix of order 976 with an
  involutive anti-automorphism ([NOTES.md](../../NOTES.md) §29).
- The symmetric quaternary form of Yamada's matrix is the new step ([YAMADA_KIYASU.md](../YAMADA_KIYASU.md#status));
  we found it in none of the sources listed there. The skew input is classical.
- Ordinary Hadamard matrices of order 7804 are classical (Yamada 1989, Theorem 4, from the same skew input); the
  symmetry is what is new.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, r, ξ, the complete certificate of the skew input, and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/7804
```

`verify/verify.py` (construction `symmetric-yamada-kiyasu`) rebuilds the skew input from its certificate and checks it
(Q + Qᵀ = 2I), computes Yamada's circulants and checks M² + N² = 1949I, builds X and Y and checks that they are
symmetric with disjoint supports covering every position, and checks H = Hᵀ, H Hᵀ = 7804 I and the SHA-256 of the
canonical CSV.

SHA-256 of `H7804.csv`: `c840b35bbe34860e218542a9ad54136d3c6aa16d5078af0e978e03f8f3062096`

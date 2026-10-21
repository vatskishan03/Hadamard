# Symmetric Hadamard matrix of order 6556

6556 = 4 · 1639, with 1639 = 11 · 149; 1639 = q + 2 with q = 1637 a prime ≡ 5 (mod 8) and (q + 3)/2 = 820.

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is the case q = 1637 of the construction in
[YAMADA_KIYASU.md](../YAMADA_KIYASU.md): Yamada's (Kiyasu's) Lemma 26 data put in quaternary form, which is symmetric
for any skew input, and Turyn's map to a real matrix.

- **The skew input.** Our skew-Hadamard matrix of order 820 ([skew/820](../../skew/820)). Normalised to first row
  all +1, its core S of order 819 has Sᵀ = −S, Se = 0 and S² = J − 819I.
- **Yamada's circulants.** M and N of order 819, with first rows β and α of Yamada's Lemma 3 over
  GF(1637²) = GF(1637)[w]/(w² − 2) and ξ = 2 + w primitive: symmetric, M of zero diagonal,
  M² + N² = 1637I.
- **The quaternary matrix.** With G = ½((M + N) + i(M − N)), P = [[1, −1], [−1, 1]] and e the all-one row of length
  1638, L = [[−1, −e, −1, ie], [−eᵀ, ½(1 + i)I ⊗ J₂ + G ⊗ P, −ieᵀ, −iI ⊗ P − S ⊗ J₂], [−1, −ie, 1, e],
  [ieᵀ, −iI ⊗ P + S ⊗ J₂, eᵀ, ½(i − 1)I ⊗ J₂ + Ḡ ⊗ P]] of order 3278 is symmetric with entries ±1, ±i and
  L L̄ = 3278I.
- **The matrix.** With L = X + iY, H = [[X + Y, X − Y], [X − Y, −X − Y]].

## Search

None: the matrix is determined by q, r, ξ and the skew input.

**Reproduction.**

| Command | Output |
|---|---|
| `python3 search/sym_yk.py certify 1639` | `symmetric/6556/certificate.json`, identical to the stored file (the program refuses to overwrite it: run it in a copy of the repository without that file) |
| `python3 verify/verify.py symmetric/6556` | `ok`, SHA-256 `e37eb01acffcd52615dcd43991f1323e91ce2fc303d9783cf3184b45a25684d5` |

## Status

**New** as a symmetric Hadamard matrix.

- The Đoković–Seberry census (2025-05-28) lists 1639 among its 595 undecided v for symmetric matrices.
  Seberry–Yamada's Corollary 27(i) in the propus array would need a skew-Hadamard matrix of order 820 with an
  involutive anti-automorphism ([NOTES.md](../../NOTES.md) §29).
- The symmetric quaternary form of Yamada's matrix is the new step ([YAMADA_KIYASU.md](../YAMADA_KIYASU.md#status));
  we found it in none of the sources listed there. The skew input is ours (2026).
- Hadamard matrices of order 6556 were already known (n = 1639 is not among the open n of Cati–Pasechnik's Table 4);
  the symmetry is what is new.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, r, ξ, the complete certificate of the skew input, and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/6556
```

`verify/verify.py` (construction `symmetric-yamada-kiyasu`) rebuilds the skew input from its certificate and checks it
(Q + Qᵀ = 2I), computes Yamada's circulants and checks M² + N² = 1637I, builds X and Y and checks that they are
symmetric with disjoint supports covering every position, and checks H = Hᵀ, H Hᵀ = 6556 I and the SHA-256 of the
canonical CSV.

SHA-256 of `H6556.csv`: `e37eb01acffcd52615dcd43991f1323e91ce2fc303d9783cf3184b45a25684d5`

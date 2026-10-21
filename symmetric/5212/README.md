# Symmetric Hadamard matrix of order 5212

5212 = 4 · 1303, with 1303 prime; 1303 = q + 2 with q = 1301 a prime ≡ 5 (mod 8) and (q + 3)/2 = 652.

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is the case q = 1301 of the construction in
[YAMADA_KIYASU.md](../YAMADA_KIYASU.md): Yamada's (Kiyasu's) Lemma 26 data put in quaternary form, which is symmetric
for any skew input, and Turyn's map to a real matrix.

- **The skew input.** Đoković's skew supplementary difference sets over ℤ₁₆₃ (1992; as stored in SageMath), a skew
  Goethals–Seidel matrix of order 652. Normalised to first row all +1, its core S of order 651 has Sᵀ = −S, Se = 0 and
  S² = J − 651I.
- **Yamada's circulants.** M and N of order 651, with first rows β and α of Yamada's Lemma 3 over
  GF(1301²) = GF(1301)[w]/(w² − 2) and ξ = 3 + w primitive: symmetric, M of zero diagonal,
  M² + N² = 1301I.
- **The quaternary matrix.** With G = ½((M + N) + i(M − N)), P = [[1, −1], [−1, 1]] and e the all-one row of length
  1302, L = [[−1, −e, −1, ie], [−eᵀ, ½(1 + i)I ⊗ J₂ + G ⊗ P, −ieᵀ, −iI ⊗ P − S ⊗ J₂], [−1, −ie, 1, e],
  [ieᵀ, −iI ⊗ P + S ⊗ J₂, eᵀ, ½(i − 1)I ⊗ J₂ + Ḡ ⊗ P]] of order 2606 is symmetric with entries ±1, ±i and
  L L̄ = 2606I.
- **The matrix.** With L = X + iY, H = [[X + Y, X − Y], [X − Y, −X − Y]].

## Search

None: the matrix is determined by q, r, ξ and the skew input.

**Reproduction.**

| Command | Output |
|---|---|
| `python3 search/sym_yk.py certify 1303` | `symmetric/5212/certificate.json`, identical to the stored file (the program refuses to overwrite it: run it in a copy of the repository without that file) |
| `python3 verify/verify.py symmetric/5212` | `ok`, SHA-256 `b7d360982d0a50243a2ac31ffd1a6f25d34567892f741c566db4011f81e33f62` |

## Status

**New** as a symmetric Hadamard matrix.

- The Đoković–Seberry census (2025-05-28) lists 1303 among its 595 undecided v for symmetric matrices.
  Seberry–Yamada's Corollary 27(i) in the propus array would need a skew-Hadamard matrix of order 652 with an
  involutive anti-automorphism ([NOTES.md](../../NOTES.md) §29).
- The symmetric quaternary form of Yamada's matrix is the new step ([YAMADA_KIYASU.md](../YAMADA_KIYASU.md#status));
  we found it in none of the sources listed there. The skew input is classical.
- Ordinary Hadamard matrices of order 5212 are classical (Yamada 1989, Theorem 4, from the same skew input); the
  symmetry is what is new.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, r, ξ, the complete certificate of the skew input, and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/5212
```

`verify/verify.py` (construction `symmetric-yamada-kiyasu`) rebuilds the skew input from its certificate and checks it
(Q + Qᵀ = 2I), computes Yamada's circulants and checks M² + N² = 1301I, builds X and Y and checks that they are
symmetric with disjoint supports covering every position, and checks H = Hᵀ, H Hᵀ = 5212 I and the SHA-256 of the
canonical CSV.

SHA-256 of `H5212.csv`: `b7d360982d0a50243a2ac31ffd1a6f25d34567892f741c566db4011f81e33f62`

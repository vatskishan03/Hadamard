# Symmetric Hadamard matrix of order 9340

9340 = 4 · 2335, with 2335 = 5 · 467; 2335 = q + 2 with q = 2333 a prime ≡ 5 (mod 8) and (q + 3)/2 = 1168.

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is the case q = 2333 of the construction in
[YAMADA_KIYASU.md](../YAMADA_KIYASU.md): Yamada's (Kiyasu's) Lemma 26 data put in quaternary form, which is symmetric
for any skew input, and Turyn's map to a real matrix.

- **The skew input.** A skew Goethals–Seidel matrix of order 292 over ℤ₇₃ found by `search/skew_exhaustive.py 73 2`
  (blocks invariant under a multiplier subgroup of order 9; skew-Hadamard matrices of order 292 are classical, e.g.
  Đoković's supplementary difference sets of 2023), doubled twice. Normalised to first row all +1, its core S of order
  1167 has Sᵀ = −S, Se = 0 and S² = J − 1167I.
- **Yamada's circulants.** M and N of order 1167, with first rows β and α of Yamada's Lemma 3 over
  GF(2333²) = GF(2333)[w]/(w² − 2) and ξ = 2 + w primitive: symmetric, M of zero diagonal,
  M² + N² = 2333I.
- **The quaternary matrix.** With G = ½((M + N) + i(M − N)), P = [[1, −1], [−1, 1]] and e the all-one row of length
  2334, L = [[−1, −e, −1, ie], [−eᵀ, ½(1 + i)I ⊗ J₂ + G ⊗ P, −ieᵀ, −iI ⊗ P − S ⊗ J₂], [−1, −ie, 1, e],
  [ieᵀ, −iI ⊗ P + S ⊗ J₂, eᵀ, ½(i − 1)I ⊗ J₂ + Ḡ ⊗ P]] of order 4670 is symmetric with entries ±1, ±i and
  L L̄ = 4670I.
- **The matrix.** With L = X + iY, H = [[X + Y, X − Y], [X − Y, −X − Y]].

## Search

The skew input comes from `search/skew_exhaustive.py 73 2` (see Construction); the matrix is determined by q, r, ξ and
that input.

**Reproduction.**

| Command | Output |
|---|---|
| `python3 search/skew_exhaustive.py 73 2` | prints skew Goethals–Seidel matrices of order 292 over ℤ₇₃ invariant under the multiplier 2 (row-sum pattern (1, 1, 17) of the blocks B, C, D), each verified as H Hᵀ = 292 I, H + Hᵀ = 2I; the input stored in `search/data/yk_skew_inputs.json` is one of the solutions, with row sums (1; 1, 1, 17) |
| `python3 search/sym_yk.py certify 2335` | `symmetric/9340/certificate.json`, identical to the stored file (the program refuses to overwrite it: run it in a copy of the repository without that file) |
| `python3 verify/verify.py symmetric/9340` | `ok`, SHA-256 `02f79d7f90c5db3fb87d0aee053a1705638729dff0dc262e8562385f70dc5ba2` |

## Status

**New** as a symmetric Hadamard matrix.

- The Đoković–Seberry census (2025-05-28) lists 2335 among its 595 undecided v for symmetric matrices.
  Seberry–Yamada's Corollary 27(i) in the propus array would need a skew-Hadamard matrix of order 1168 with an
  involutive anti-automorphism ([NOTES.md](../../NOTES.md) §29).
- The symmetric quaternary form of Yamada's matrix is the new step ([YAMADA_KIYASU.md](../YAMADA_KIYASU.md#status));
  we found it in none of the sources listed there. The skew input is classical.
- Ordinary Hadamard matrices of order 9340 are classical (Yamada 1989, Theorem 4, from the same skew input); the
  symmetry is what is new.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, r, ξ, the complete certificate of the skew input, and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/9340
```

`verify/verify.py` (construction `symmetric-yamada-kiyasu`) rebuilds the skew input from its certificate and checks it
(Q + Qᵀ = 2I), computes Yamada's circulants and checks M² + N² = 2333I, builds X and Y and checks that they are
symmetric with disjoint supports covering every position, and checks H = Hᵀ, H Hᵀ = 9340 I and the SHA-256 of the
canonical CSV.

SHA-256 of `H9340.csv`: `02f79d7f90c5db3fb87d0aee053a1705638729dff0dc262e8562385f70dc5ba2`

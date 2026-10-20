# Symmetric Hadamard matrix of order 4732

4732 = 4 · 1183, with 1183 = 7 · 13²; 1183 = q + 2 with q = 1181 a prime ≡ 5 (mod 8) and (q + 3)/2 = 592.

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is the case q = 1181 of the construction in
[YAMADA_KIYASU.md](../YAMADA_KIYASU.md): Yamada's (Kiyasu's) Lemma 26 data put in quaternary form, which is symmetric for
any skew input, and Turyn's map to a real matrix.

- **The skew input.** A skew Goethals–Seidel matrix of order 148 over ℤ₃₇ found by `search/skew_exhaustive.py 37 10` (blocks invariant under a multiplier subgroup of order 3; skew-Hadamard matrices of order 148 are classical, e.g. Đoković's supplementary difference sets of 1992), doubled twice. Normalised to first row all +1, its core S of order 591 has
  Sᵀ = −S, Se = 0 and S² = J − 591I.
- **Yamada's circulants.** M and N of order 591, with first rows β and α of Yamada's Lemma 3 over
  GF(1181²) = GF(1181)[w]/(w² − 2) and ξ = 3 + w primitive: symmetric, M of zero diagonal,
  M² + N² = 1181I.
- **The quaternary matrix.** With G = ½((M + N) + i(M − N)), P = [[1, −1], [−1, 1]] and e the all-one row of length 1182,
  L = [[−1, −e, −1, ie], [−eᵀ, ½(1 + i)I ⊗ J₂ + G ⊗ P, −ieᵀ, −iI ⊗ P − S ⊗ J₂], [−1, −ie, 1, e],
  [ieᵀ, −iI ⊗ P + S ⊗ J₂, eᵀ, ½(i − 1)I ⊗ J₂ + Ḡ ⊗ P]] of order 2366 is symmetric with entries ±1, ±i and
  L L̄ = 2366I.
- **The matrix.** With L = X + iY, H = [[X + Y, X − Y], [X − Y, −X − Y]].

## Search

None: the matrix is determined by q, r, ξ and the skew input. `python3 search/sym_yk.py certify 1183` writes the
certificate.

## Status

**New** as a symmetric Hadamard matrix.

- The Đoković–Seberry census (2025-05-28) lists 1183 among its 595 undecided v for symmetric matrices. Seberry–Yamada's
  Corollary 27(i) in the propus array would need a symmetric conference matrix of order 592 or a skew-Hadamard
  matrix of order 592 with an involutive anti-automorphism ([NOTES.md](../../NOTES.md) §29).
- The symmetric quaternary form of Yamada's matrix is the new step ([YAMADA_KIYASU.md](../YAMADA_KIYASU.md#status)); we
  found it in none of the sources listed there. The skew input is classical.
- Ordinary Hadamard matrices of order 4732 are classical (Yamada 1989, Theorem 4, from the same skew input); the
  symmetry is what is new.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, r, ξ, the complete certificate of the skew input, and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/4732
```

`verify/verify.py` (construction `symmetric-yamada-kiyasu`) rebuilds the skew input from its certificate and checks it
(Q + Qᵀ = 2I), computes Yamada's circulants and checks M² + N² = 1181I, builds X and Y and checks that they are symmetric
with disjoint supports covering every position, and checks H = Hᵀ, H Hᵀ = 4732 I and the SHA-256 of the canonical CSV.

SHA-256 of `H4732.csv`: `c1926f9f339dff11414843629092b05a378025cdd54aa01226fe86a522893858`

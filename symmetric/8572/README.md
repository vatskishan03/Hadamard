# Symmetric Hadamard matrix of order 8572

8572 = 4 · 2143, with 2143 prime; 2143 = q + 2 with q = 2141 a prime ≡ 5 (mod 8) and (q + 3)/2 = 1072.

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is the case q = 2141 of the construction in
[YAMADA_KIYASU.md](../YAMADA_KIYASU.md): Yamada's (Kiyasu's) Lemma 26 data put in quaternary form, which is symmetric for
any skew input, and Turyn's map to a real matrix.

- **The skew input.** Đoković's skew supplementary difference sets over ℤ₆₇ (1992; as stored in SageMath), a skew Goethals–Seidel matrix of order 268, doubled twice. Normalised to first row all +1, its core S of order 1071 has
  Sᵀ = −S, Se = 0 and S² = J − 1071I.
- **Yamada's circulants.** M and N of order 1071, with first rows β and α of Yamada's Lemma 3 over
  GF(2141²) = GF(2141)[w]/(w² − 2) and ξ = 2 + w primitive: symmetric, M of zero diagonal,
  M² + N² = 2141I.
- **The quaternary matrix.** With G = ½((M + N) + i(M − N)), P = [[1, −1], [−1, 1]] and e the all-one row of length 2142,
  L = [[−1, −e, −1, ie], [−eᵀ, ½(1 + i)I ⊗ J₂ + G ⊗ P, −ieᵀ, −iI ⊗ P − S ⊗ J₂], [−1, −ie, 1, e],
  [ieᵀ, −iI ⊗ P + S ⊗ J₂, eᵀ, ½(i − 1)I ⊗ J₂ + Ḡ ⊗ P]] of order 4286 is symmetric with entries ±1, ±i and
  L L̄ = 4286I.
- **The matrix.** With L = X + iY, H = [[X + Y, X − Y], [X − Y, −X − Y]].

## Search

None: the matrix is determined by q, r, ξ and the skew input. `python3 search/sym_yk.py certify 2143` writes the
certificate.

## Status

**New** as a symmetric Hadamard matrix.

- The Đoković–Seberry census (2025-05-28) lists 2143 among its 595 undecided v for symmetric matrices. Seberry–Yamada's
  Corollary 27(i) in the propus array would need a symmetric conference matrix of order 1072 or a skew-Hadamard
  matrix of order 1072 with an involutive anti-automorphism ([NOTES.md](../../NOTES.md) §29).
- The symmetric quaternary form of Yamada's matrix is the new step ([YAMADA_KIYASU.md](../YAMADA_KIYASU.md#status)); we
  found it in none of the sources listed there. The skew input is classical.
- Ordinary Hadamard matrices of order 8572 are classical (Yamada 1989, Theorem 4, from the same skew input); the
  symmetry is what is new.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, r, ξ, the complete certificate of the skew input, and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/8572
```

`verify/verify.py` (construction `symmetric-yamada-kiyasu`) rebuilds the skew input from its certificate and checks it
(Q + Qᵀ = 2I), computes Yamada's circulants and checks M² + N² = 2141I, builds X and Y and checks that they are symmetric
with disjoint supports covering every position, and checks H = Hᵀ, H Hᵀ = 8572 I and the SHA-256 of the canonical CSV.

SHA-256 of `H8572.csv`: `bdacc1a59d66fb62de4d1812c3bcb1402233de22869ecd375fee58498bbc5cf1`

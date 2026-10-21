# Symmetric Hadamard matrix of order 2812

2812 = 4 · 703, with 703 = 19 · 37; 703 = q + 2 with q = 701 a prime ≡ 5 (mod 8) and (q + 3)/2 = 352.

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is the case q = 701 of the construction in
[YAMADA_KIYASU.md](../YAMADA_KIYASU.md): Yamada's (Kiyasu's) Lemma 26 data put in quaternary form, which is symmetric
for any skew input, and Turyn's map to a real matrix.

- **The skew input.** Paley's skew-Hadamard matrix I + C of order 44 (p = 43), doubled three times
  ([[S, S], [−Sᵀ, Sᵀ]]). Normalised to first row all +1, its core S of order 351 has Sᵀ = −S,
  Se = 0 and S² = J − 351I.
- **Yamada's circulants.** M and N of order 351, with first rows β and α of Yamada's Lemma 3 over
  GF(701²) = GF(701)[w]/(w² − 2) and ξ = 2 + w primitive: symmetric, M of zero diagonal,
  M² + N² = 701I.
- **The quaternary matrix.** With G = ½((M + N) + i(M − N)), P = [[1, −1], [−1, 1]] and e the all-one row of
  length 702, L = [[−1, −e, −1, ie], [−eᵀ, ½(1 + i)I ⊗ J₂ + G ⊗ P, −ieᵀ, −iI ⊗ P − S ⊗ J₂], [−1, −ie, 1, e],
  [ieᵀ, −iI ⊗ P + S ⊗ J₂, eᵀ, ½(i − 1)I ⊗ J₂ + Ḡ ⊗ P]] of order 1406 is symmetric with entries ±1, ±i and
  L L̄ = 1406I.
- **The matrix.** With L = X + iY, H = [[X + Y, X − Y], [X − Y, −X − Y]].

## Search

None: the matrix is determined by q, r, ξ and the skew input. `python3 search/sym_yk.py certify 703` writes the
certificate.

## Status

**Known** as a symmetric Hadamard matrix.

- Spence (An infinite family of Williamson matrices, J. Austral. Math. Soc. Ser. A 24 (1977) 252–256),
  Corollary 3 with q = 37 and r = 1, gives a symmetric Hadamard matrix of order 2qʳ(q + 1) = 2812. The matrix
  here is a different construction, Yamada's matrix in symmetric quaternary form
  ([YAMADA_KIYASU.md](../YAMADA_KIYASU.md)).
- The Đoković–Seberry census (2025-05-28) lists v = 703 as undecided; its series omit Spence's corollary.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, r, ξ, the complete certificate of the skew input, and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/2812
```

`verify/verify.py` (construction `symmetric-yamada-kiyasu`) rebuilds the skew input from its certificate and checks it
(Q + Qᵀ = 2I), computes Yamada's circulants and checks M² + N² = 701I, builds X and Y and checks that they are
symmetric with disjoint supports covering every position, and checks H = Hᵀ, H Hᵀ = 2812 I and the SHA-256 of the
canonical CSV.

SHA-256 of `H2812.csv`: `490f0eec754029801448f10ddddee51de4971e227972fd30983988d8a306688e`

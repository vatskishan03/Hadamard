# Symmetric Hadamard matrices of order 2n(n − 1) from a symmetric conference matrix of order n

**Theorem.** If a symmetric conference matrix of order n exists, then a symmetric Hadamard matrix of order 2n(n − 1)
exists.

With Paley's conference matrices (n = q + 1, q ≡ 1 mod 4 a prime power) this gives symmetric Hadamard matrices of order
2q(q + 1) = 4v, v = q(q + 1)/2. Five of these v are undecided in the Đoković–Seberry census (2025-05-28): 153, 325, 435,
1225 and 1891 (q = 17, 25, 29, 49, 61), the orders [612](612), [1300](1300), [1740](1740), [4900](4900) and
[7564](7564). The construction is the symmetric form of a complex Hadamard matrix of Kharaghani and Seberry, so its
ingredients are classical ([Status](#status)).

The construction was found in round 13 of the exploration of 18–19 October 2026 (the recursive-constructions lens), and
the proof below was checked by hand. Its builder in [verify/verify.py](../verify/verify.py) (construction
`symmetric-conference-plug`) is an implementation independent of the agent's.

## Proof

Let C be a symmetric conference matrix of order n: C = Cᵀ, zero diagonal, ±1 elsewhere, C² = (n − 1)I. Put m = n − 1.
Multiplying rows and columns by the signs of the first row gives an equivalent symmetric conference matrix
C = [[0, eᵀ], [e, Q]] with e the all-one vector. Then (C²)₀ⱼ = 0 gives Q𝟙 = 0, so JQ = QJ = 0, and the lower block of
C² gives Q² = mI − J.

Let J = J_m and

    X = [[J, J], [J, −J]],   Y = [[Q + I, Q − I], [Q − I, −Q − I]].

Both are symmetric of order 2m with entries ±1 (Q has zero diagonal). Using J² = mJ and JQ = QJ = 0:

- X² = I₂ ⊗ 2mJ;
- Y² = I₂ ⊗ (2(m + 1)I − 2J): the diagonal blocks are (Q + I)² + (Q − I)² = 2Q² + 2I, and the off-diagonal blocks
  (Q + I)(Q − I) − (Q − I)(Q + I) vanish;
- XY = [[0, −2J], [2J, 0]] = −YX.

Let

    H = I_n ⊗ X + C ⊗ Y.

The diagonal blocks of H are X and the others are C_ij Y with C_ij = ±1, so H has entries ±1. It is symmetric because
X, Y and C are. Finally

    H² = I ⊗ X² + C ⊗ (XY + YX) + C² ⊗ Y² = I ⊗ (X² + mY²) = I ⊗ I₂ ⊗ (2mJ + 2m(m + 1)I − 2mJ) = 2n(n − 1) I.

**Complex form.** G = I_n ⊗ J_m + C ⊗ (Q + iI) has entries ±1, ±i, and Gᵀ = G. Its real and imaginary parts
M = I ⊗ J + C ⊗ Q and E = C ⊗ I commute and have disjoint supports, and GG* = n(n − 1)I. H is, up to a permutation of
rows and columns, Turyn's symmetric image [[M + E, M − E], [M − E, −M − E]] of G (Seberry–Yamada 1992, Lemma 11.3 with
a = b = 1).

**The plug pair.** In the plug-ins of [README](README.md) (4828 and the m(m − 1), m(m + 3) constructions), the design
(P, N) anticommutes and the plugs commute. Here the design (I, C) commutes and the plugs X, Y anticommute. Within this
scheme the plugs cannot be smaller (the round-13 agent's bound, re-derived here). Take an outer conference matrix of
order n and symmetric ±1 plugs of order 2k with XY + YX = 0 and X² + (n − 1)Y² = n · 2k · I, and let the order 2nk
be 4v with v odd, so that k is odd.
- The off-diagonal entries of X² and Y² are even, of modulus at most 2k, and (X²)ᵢⱼ = −(n − 1)(Y²)ᵢⱼ.
- If they all vanish, X and Y are symmetric Hadamard matrices of order 2k ≡ 2 (mod 4). Then k = 1, which is Paley II
  (v = n/2).
- Otherwise 2(n − 1) ≤ 2k.
So apart from Paley II the scheme reaches only v ≥ n(n − 1)/2, and k = n − 1 is the case built here.

## Reach

Paley's conference matrices give v = q(q + 1)/2 for prime powers q ≡ 1 (mod 4) (`python3 search/sym_conference_plug.py
targets`):

| q | n | v | order | census | here |
|---:|---:|---:|---:|---|---|
| 5, 9, 13 | 6, 10, 14 | 15, 45, 91 | 60, 180, 364 | known | |
| 17 | 18 | 153 = 3² · 17 | [612](612) | undecided | New, but obtainable classically |
| 25 | 26 | 325 = 5² · 13 | [1300](1300) | undecided | New, but obtainable classically |
| 29 | 30 | 435 = 3 · 5 · 29 | [1740](1740) | undecided | New, but obtainable classically |
| 37 | 38 | 703 = 19 · 37 | 2812 | undecided | already [2812](2812) (Yamada–Kiyasu) |
| 41 | 42 | 861 | 3444 | known | |
| 49 | 50 | 1225 = 5² · 7² | [4900](4900) | undecided | New, but obtainable classically |
| 53 | 54 | 1431 | 5724 | known | |
| 61 | 62 | 1891 = 31 · 61 | [7564](7564) | undecided | New, but obtainable classically |

The next Paley value, q = 73, gives v = 2701, beyond the census. The other known symmetric conference matrices give
nothing new below 2500: Mathon's order 46 gives v = 1035 (known), and Verma's order 86 gives 3655. Kharaghani and
Seberry's recursion G_{j+1} = I ⊗ (J ⊗ B_j) + C ⊗ (Q ⊗ B_j + I ⊗ A_j) gives symmetric orders 2n(n − 1)^j as well. For
v < 2500 these are v = 75, 375, 405, 1183 and 1875, which are known or already here ([4732](4732)).

## Status

The five orders are **New, but obtainable classically**. All five v are undecided in the census, and we found no
symmetric Hadamard matrix of these orders and no statement of the construction (audit of 19 October 2026,
[NOVELTY.md](../NOVELTY.md)). Every ingredient was published before 2026:
- Paley's conference matrices;
- the complex Hadamard matrix G of order n(n − 1). It is Kharaghani and Seberry's Theorem 11.14 with the skew-type
  complex Hadamard matrix I + iC, A = J and B = I − iQ, which satisfy AB* = BA* and AA* + mBB* = m(m + 1)I
  (Seberry–Yamada 1992, Theorem 11.14 and Corollary 11.15). Corollary 11.20 and Table 11.1 there list the orders
  w(w − 1) for w = 18, 26, 30, 38, 50, 62;
- Turyn's map, which is symmetric for a complex-symmetric input (Seberry–Yamada 1992, Lemma 11.3).

Only the observation that G is complex-symmetric, with commuting real and imaginary parts, is new to the sources we
read. Ordinary Hadamard matrices of these orders were known, for instance from Whiteman's Williamson matrices of order
q(q + 1)/2 (1976), which Seberry (J. Austral. Math. Soc. A 29 (1980), final remarks) mentions. The symmetry is what is
new.

## Reproduction and verification

```
python3 search/sym_conference_plug.py targets                 # v = q(q + 1)/2 in the census
python3 search/sym_conference_plug.py certify 17 25 29 49 61  # writes the five certificates
python3 verify/verify.py symmetric/612                        # rebuilds H and checks it
```

`verify/verify.py` (construction `symmetric-conference-plug`) builds Paley's conference matrix over GF(q). For q = 25
and 49 the field is GF(p)[x]/(poly), with the primitive polynomial given in the certificate. The builder checks that
C is a normalised symmetric conference matrix, that Q𝟙 = 0 and XY + YX = 0. It then assembles H and checks H = Hᵀ,
HHᵀ = nI and the SHA-256 of the canonical CSV.

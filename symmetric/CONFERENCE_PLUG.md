# Symmetric Hadamard matrices of order 2n(n − 1) from a symmetric conference matrix of order n

**Theorem** (Spence 1977, Corollary 2 with r = 1). If a symmetric conference matrix of order n exists, then a
symmetric Hadamard matrix of order 2n(n − 1) exists.

Spence (An infinite family of Williamson matrices, J. Austral. Math. Soc. Ser. A 24 (1977) 252–256) proved more: a
conference matrix of order n + 1 with n ≡ 1 (mod 4) gives a symmetric Hadamard matrix of order 2nʳ(n + 1) for every
r ≥ 0 (Corollary 2), in particular for every prime power n = q ≡ 1 (mod 4) (Corollary 3). The construction below,
H = I ⊗ X + C ⊗ Y, was found independently; it gives the case r = 1 in another form, and whether its matrices are
equivalent to Spence's was not checked.

With Paley's conference matrices (n = q + 1, q ≡ 1 mod 4 a prime power) it gives symmetric Hadamard matrices of order
2q(q + 1) = 4v, v = q(q + 1)/2. Five of these v are undecided in the Đoković–Seberry census (2025-05-28), whose series
omit Spence's corollary: 153, 325, 435, 1225 and 1891 (q = 17, 25, 29, 49, 61), the orders [612](612), [1300](1300),
[1740](1740), [4900](4900) and [7564](7564). They are Known ([Status](#status)). The builder in
[verify/verify.py](../verify/verify.py) (construction `symmetric-conference-plug`) rebuilds and checks them.

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
scheme the plugs cannot be smaller. Take an outer conference matrix of order n and symmetric ±1 plugs of order 2k with
XY + YX = 0 and X² + (n − 1)Y² = n · 2k · I, and let the order 2nk be 4v with v odd, so that k is odd.
- The off-diagonal entries of X² and Y² are even, of modulus at most 2k, and (X²)ᵢⱼ = −(n − 1)(Y²)ᵢⱼ.
- If they all vanish, X and Y are symmetric Hadamard matrices of order 2k ≡ 2 (mod 4). Then k = 1, which is Paley II
  (v = n/2).
- Otherwise 2(n − 1) ≤ 2k.
So apart from Paley II the scheme reaches only v ≥ n(n − 1)/2, and k = n − 1 is the case built here.

## Reach

Paley's conference matrices give v = q(q + 1)/2 for prime powers q ≡ 1 (mod 4)
(`python3 search/sym_conference_plug.py targets`):

| q | n | v | order | census | here |
|---:|---:|---:|---:|---|---|
| 5, 9, 13 | 6, 10, 14 | 15, 45, 91 | 60, 180, 364 | known | |
| 17 | 18 | 153 = 3² · 17 | [612](612) | undecided | Known (Spence 1977) |
| 25 | 26 | 325 = 5² · 13 | [1300](1300) | undecided | Known (Spence 1977) |
| 29 | 30 | 435 = 3 · 5 · 29 | [1740](1740) | undecided | Known (Spence 1977) |
| 37 | 38 | 703 = 19 · 37 | [2812](2812) | undecided | Known (Spence 1977); here from Yamada–Kiyasu |
| 41 | 42 | 861 | 3444 | known | |
| 49 | 50 | 1225 = 5² · 7² | [4900](4900) | undecided | Known (Spence 1977) |
| 53 | 54 | 1431 | 5724 | known | |
| 61 | 62 | 1891 = 31 · 61 | [7564](7564) | undecided | Known (Spence 1977) |

The next Paley value, q = 73, gives v = 2701, beyond the census. The other known symmetric conference matrices give
nothing new below 2500: Mathon's order 46 gives v = 1035 (known), Gritsenko's order 66 (2021) gives v = 2145 (known),
and Verma's order 86 gives 3655. Spence's recursion (Corollary 2) gives orders 2n(n − 1)^j, and the builder implements
only the case j = 1. For v < 2500 the cases with j ≥ 2 are v = 75, 375, 405, 1183 and 1875, which are known or already
here (v = 1183 is [4732](4732), j = 2, built by Yamada–Kiyasu).

## Status

The five orders are **Known**: Spence's Corollary 3 with r = 1 and q = 17, 25, 29, 49 and 61 gives symmetric Hadamard
matrices of these orders. The census lists the five v as undecided because its series omit the corollary. The
matrices here were built with the construction above, found independently. The corollary also covers v = 703
(q = 37, r = 1) and v = 1183 (q = 13, r = 2), the orders [2812](2812) and [4732](4732), which are here from
Yamada–Kiyasu.

The construction is also the symmetric form of a complex Hadamard matrix of Kharaghani and Seberry (Seberry–Yamada
1992, Theorem 11.14 and Corollary 11.15) under Turyn's map (Lemma 11.3 there). Ordinary Hadamard matrices of these
orders were known before Spence's paper, from the Williamson matrices of order q(q + 1)/2 of Wallis (1973) and
Whiteman (1976).

The sources searched are listed in [NOVELTY.md](../NOVELTY.md).

## Reproduction and verification

```
python3 search/sym_conference_plug.py targets                 # v = q(q + 1)/2 in the census
python3 search/sym_conference_plug.py certify 17 25 29 49 61  # writes the five certificates
python3 verify/verify.py symmetric/612                        # rebuilds H and checks it
```

`certify` writes a certificate only into a directory that has none: if `symmetric/<order>/certificate.json` exists it
stops with `<directory> exists`. The five certificates are committed, so to regenerate one, move its
`certificate.json` away first and compare the new file with the old.

`verify/verify.py` (construction `symmetric-conference-plug`) builds Paley's conference matrix over GF(q). For q = 25
and 49 the field is GF(p)[x]/(poly), with the primitive polynomial given in the certificate. The builder checks that C
is a normalised symmetric conference matrix, that Q𝟙 = 0 and XY + YX = 0. It then assembles H and checks H = Hᵀ,
H Hᵀ = 2n(n − 1)I and the SHA-256 of the canonical CSV.

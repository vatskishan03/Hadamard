# Symmetric Hadamard matrices of order 4(q + 2) from any skew-Hadamard matrix of order (q + 3)/2

**Theorem.** Let q ≡ 5 (mod 8) be a prime power. If a skew-Hadamard matrix of order (q + 3)/2 exists, then a symmetric
Hadamard matrix of order 4(q + 2) exists.

Yamada (J. Austral. Math. Soc. Ser. A 46 (1989), Theorem 4, credited to Kiyasu) built ordinary Hadamard matrices of
order 4(q + 2) from the same input ([search/NOTES.md](../search/NOTES.md) §13, [hadamard/15292](../hadamard/15292)). The
symmetric route of the literature, Seberry–Yamada's Corollary 27(i) in Seberry–Balonin's propus array, needs a
symmetric conference core, or (NOTES.md §29) a skew core with an involutive anti-automorphism, so that the
Williamson-type matrices have a repeated block. Here the repeated block is not needed: Yamada's Lemma 26 data are put in
quaternary form, as a symmetric matrix L with entries ±1, ±i and L L̄ = 2(q + 2)I, and Turyn's map
[[X + Y, X − Y], [X − Y, −X − Y]] of L = X + iY is a symmetric Hadamard matrix (the map is classical: Seberry–Yamada
1992, Lemma 11.3 with a = b = 1, symmetric whenever L is).

The theorem settles 27 values of the Đoković–Seberry census (2025-05-28), from v = 295 to 2335. For 295 the skew-propus
route of NOTES.md §29 had failed (no skew-propus quadruple of length 37, search/NOTES.md §49).

| Order | v = q + 2 | (q + 3)/2 | Skew-Hadamard input | Status |
|------:|---:|---:|---|---|
| [1180](1180) | 295 = 5 · 59 | 148 | Goethals–Seidel over ℤ₃₇ (our search; Đoković 1992) | New |
| [1276](1276) | 319 = 11 · 29 | 160 | Paley (p = 19), doubled three times | New |
| [1564](1564) | 391 = 17 · 23 | 196 | Goethals–Seidel over ℤ₄₉ (our search; Đoković 1992) | New |
| [1852](1852) | 463 | 232 | good matrices of order 29 (Szekeres 1988), doubled | New |
| [2620](2620) | 655 = 5 · 131 | 328 | Paley (p = 163), doubled | New |
| [2716](2716) | 679 = 7 · 97 | 340 | Whiteman's complementary difference sets in GF(13²) (1971) | New |
| [2812](2812) | 703 = 19 · 37 | 352 | Paley (p = 43), doubled three times | New |
| [3292](3292) | 823 | 412 | Goethals–Seidel over ℤ₁₀₃ (Đoković 1994) | New |
| [3772](3772) | 943 = 23 · 41 | 472 | Goethals–Seidel over ℤ₅₉ (Fletcher–Koukouvinos–Seberry 2004), doubled | New |
| [4252](4252) | 1063 | 532 | Goethals–Seidel over ℤ₁₃₃ (our search; Đoković 1992) | New |
| [4732](4732) | 1183 = 7 · 13² | 592 | Goethals–Seidel over ℤ₃₇, doubled twice | New |
| [4924](4924) | 1231 | 616 | Paley (p = 307), doubled | New |
| [5212](5212) | 1303 | 652 | Goethals–Seidel over ℤ₁₆₃ (Đoković 1992) | New |
| [5980](5980) | 1495 = 5 · 13 · 23 | 748 | Szekeres's complementary difference sets in GF(373) (1969) | New |
| [6396](6396) | 1599 = 3 · 13 · 41 | 800 | Paley (p = 199), doubled twice | New |
| [6556](6556) | 1639 = 11 · 149 | 820 | our skew matrix [skew/820](../skew/820) | New |
| [6684](6684) | 1671 = 3 · 557 | 836 | our skew matrix [skew/836](../skew/836) | New |
| [7804](7804) | 1951 | 976 | Paley (p = 487), doubled | New |
| [7900](7900) | 1975 = 5² · 79 | 988 | Goethals–Seidel over ℤ₂₄₇ (Đoković 2008) | New |
| [7996](7996) | 1999 | 1000 | Paley (p = 499), doubled | New |
| [8284](8284) | 2071 = 19 · 109 | 1036 | our skew matrix [skew/1036](../skew/1036) | New |
| [8572](8572) | 2143 | 1072 | Goethals–Seidel over ℤ₆₇ (Đoković 1992), doubled twice | New |
| [8892](8892) | 2223 = 3² · 13 · 19 | 1112 | Szekeres's complementary difference sets in GF(277), doubled | New |
| [8956](8956) | 2239 | 1120 | Paley (p = 139), doubled three times | New |
| [9084](9084) | 2271 = 3 · 757 | 1136 | Paley (p = 283), doubled twice | New |
| [9180](9180) | 2295 = 3³ · 5 · 17 | 1148 | our skew matrix [skew/1148](../skew/1148) | New |
| [9340](9340) | 2335 = 5 · 467 | 1168 | Goethals–Seidel over ℤ₇₃, doubled twice | New |

## Proof

Let n = (q + 1)/2 (≡ 3 mod 4) and v = q + 2 = 2n + 1.

- **The skew core.** Normalise the skew-Hadamard matrix of order n + 1 so that its first row is all +1; its core S
  (delete the first row and column, subtract I) satisfies Sᵀ = −S, Se = 0 and S² = −SSᵀ = J − nI.
- **Yamada's circulants.** M and N, symmetric circulants of order n with first rows β and α of Yamada's Lemma 3 (from a
  primitive element of GF(q²), as in [hadamard/15292](../hadamard/15292)): M has zero diagonal and entries ±1 off it,
  N has entries ±1, MN = NM and M² + N² = qI.

Put G = ½((M + N) + i(M − N)), Ḡ its complex conjugate, P = [[1, −1], [−1, 1]], J₂ the all-one matrix of order 2, e
the all-one row of length 2n, and index the cores by (j, a), j < n, a ∈ {0, 1}. In blocks of orders 1, 2n, 1, 2n:

```
     [ −1      −e                         −1       i e                      ]
L =  [ −eᵀ     ½(1 + i) I ⊗ J₂ + G ⊗ P    −i eᵀ    −i I ⊗ P − S ⊗ J₂         ]
     [ −1      −i e                        1       e                        ]
     [ i eᵀ    −i I ⊗ P + S ⊗ J₂           eᵀ      ½(i − 1) I ⊗ J₂ + Ḡ ⊗ P   ]
```

**Entries and symmetry.** Off the diagonal, G_jk = ½((M + N) + i(M − N))_jk is one of ±1, ±i because M_jk, N_jk = ±1;
on it, G_jj = ½(1 − i)N_jj because M_jj = 0, and ½(1 + i)J₂ + ½(1 − i)N_jj P is [[1, i], [i, 1]] or [[i, 1], [1, i]]
(likewise for the last diagonal block). The other core blocks have entries −iP_aa′ (j = k, as S_jj = 0) or ∓S_jk.
L is symmetric: G is, and (−iI ⊗ P − S ⊗ J₂)ᵀ = −iI ⊗ P + S ⊗ J₂ because Sᵀ = −S.

**L L̄ = 2vI.** J₂P = PJ₂ = 0 and e(Y ⊗ P) = 0, so L and L̄ preserve the subspace 𝒱 spanned by the vectors x ⊗ (1, −1) on
either core and the subspace 𝒰 spanned by the borders and the vectors x ⊗ (1, 1); the two are orthogonal and fill the
space. In the basis x ⊗ (1, −1)/√2 of 𝒱 (left core, right core),

```
L|𝒱 = 2 [[G, −iI], [−iI, Ḡ]],    L|𝒱 · L̄|𝒱 = 4 [[GḠ + I, 0], [0, ḠG + I]] = 4(q/2 + 1) I = 2v I,
```

since GḠ = ¼((M + N)² + (M − N)²) = ½(M² + N²) = (q/2)I. In the basis (0_L, x ⊗ (1, 1)/√2, 0_R, x ⊗ (1, 1)/√2) of 𝒰,
with 1 the all-one column of length n,

```
        [ −1      −√2 1ᵀ     −1       √2 i 1ᵀ ]
L|𝒰 =   [ −√2 1   (1 + i)I   −√2 i 1   −2S     ]
        [ −1      −√2 i 1ᵀ    1       √2 1ᵀ   ]
        [ √2 i 1   2S         √2 1     (i − 1)I ]
```

and a block computation with Sᵀ = −S, 1ᵀS = 0 and S² = J − nI gives L|𝒰 · L̄|𝒰 = (4n + 2)I = 2vI. For example the
(core, core) blocks are 2J + (1 + i)(1 − i)I + 2J − 4S² = (4n + 2)I on the left and 2J − 4S² + 2J + (i − 1)(−i − 1)I =
(4n + 2)I on the right, the (left core, right core) block is 2iJ − 2(1 + i)S − 2iJ + 2(1 + i)S = 0, and the border rows
give 1 + 2n + 1 + 2n = 4n + 2 and, against the cores, √2 1ᵀ(1 − (1 − i) − i) = 0. L L̄ is Hermitian (L is symmetric), so
the remaining blocks follow.

**The real matrix.** Write L = X + iY with X, Y real. They are symmetric with disjoint supports covering every position,
and L L̄ = X² + Y² + i(YX − XY) = 2vI gives X² + Y² = 2vI and XY = YX. So

```
H = [[X + Y,  X − Y],
     [X − Y,  −X − Y]]
```

has entries ±1, is symmetric, and H² = [[2(X² + Y²), 2(YX − XY)], [2(XY − YX), 2(X² + Y²)]] = 4vI. ∎

Remarks.

- L is Yamada's Lemma 26 data U = (I, 0, S, −S), V = (M, N, I, I) in quaternary form. With the bordered matrices
  W_k = [[c_k, ε_k e], [ε_k eᵀ, U_k ⊗ J₂ + V_k ⊗ P]], corners c = (1, −1, 1, 1) and borders ε = (−1, −1, 1, 1)
  (W₄ = W₃ᵀ, not symmetric), X′ = ½((W₁ + W₂) + i(W₁ − W₂)) and Y′ = ½((W₃ + W₄) + i(W₃ − W₄)), the quaternary
  Hadamard matrix [[X′, Y′], [−iY′, iX̄′]] times a monomial matrix is L: swap the two border columns and the two
  columns of each right-core pair, then multiply the columns by −1, 1, i, i on the four blocks (checked for q = 13
  and 37).
- Only Sᵀ = −S, Se = 0, S² = J − nI and the pair M, N are used: no anti-automorphism of S, and any symmetric
  commuting M (zero diagonal) and N with M² + N² = qI would do. q ≡ 5 (mod 8) is forced by the existence of a
  skew-Hadamard matrix of order n + 1 = (q + 3)/2.

## Reach

Of the census's undecided v, 35 are q + 2 with q ≡ 5 (mod 8) a prime power (189, 473 and 1139, settled by Đoković, are
not of this form). One of them, 183, was settled here earlier by a search ([732](732)); the theorem reaches it too,
with a skew-Hadamard matrix of order 92. Of the other 34, a skew-Hadamard matrix of order (v + 1)/2 is known for 27:

- Paley's matrices (p ≡ 3 mod 4), doubled: v = 319, 655, 703, 1231, 1599, 1951, 1999, 2239, 2271;
- Szekeres's (1969) and Whiteman's (1971) complementary difference sets: v = 679, 1495, 2223;
- Goethals–Seidel matrices from published data: Szekeres's good matrices of order 29 (1988), Fletcher, Koukouvinos and
  Seberry's matrix of order 236 (2004) and Đoković's supplementary difference sets for 67, 103, 163 and 247 (1992–2008),
  all as stored in SageMath and checked here: v = 463, 943, 2143, 823, 1303, 1975;
- Goethals–Seidel matrices found by `search/skew_exhaustive.py` in at most a few minutes (orders 148, 196, 292 and
  532, all classical: Đoković's supplementary difference sets of 1992 and 2023 give them too): v = 295, 1183, 391,
  2335, 1063;
- our skew matrices of orders 820, 836, 1036 and 1148: v = 1639, 1671, 2071, 2295.

The other seven, v = 1711, 1879, 1903, 2215, 2359, 2383 and 2391, need skew-Hadamard matrices of orders 856 = 8 · 107,
940 = 4 · 235, 952 = 8 · 119, 1108 = 4 · 277, 1180 = 4 · 295, 1192 = 8 · 149 and 1196 = 4 · 299, all undecided in the
skew census (as are the skew orders 428, 476 and 596 whose doublings would serve). A skew matrix of any of these orders
would decide one more symmetric value.

## Status

- All 27 orders are **New**: their census values are undecided, and we found the symmetric form of Yamada's matrix
  nowhere (audit of 17 October 2026). Ordinary Hadamard matrices of all these orders are classical (Yamada 1989,
  Theorem 4, announced in Proc. Japan Acad. 63 (1987); Đoković, arXiv:1008.2043, applies it to 823, 1063, 1303, 2143
  and 2335); the symmetry is what is new.
  - The symmetric case published so far is Paley's skew core: Seberry–Yamada's Corollary 27(i) (JCMCC 7 (1990); 1992
    survey, Corollary 8.7; book of 2020, Corollary 6.3) gives Williamson-type matrices with a repeated block when
    (q + 1)/2 is a prime power, and Seberry–Balonin's propus array (2017; the order corrected by Balonin–Đoković,
    arXiv:1803.08787) makes them symmetric. For none of the 27 values is (q + 1)/2 a prime power.
  - Not read: Kiyasu's book (1980, in Japanese) and his 1974 series, which state the ordinary result; the full text
    of Seberry–Yamada (2020), chapters 3 and 6 and pp. 211–214; Miyamoto (1991). The census of 2025, co-authored by
    Seberry, lists all 27 values as undecided.
- **Zaitsev et al.** The page [Quaternary Hadamard matrices](https://www.cs.uleth.ca/~zaitsevv/quaternary.html) of
  Zaitsev, Kharaghani, Mohammadian and Tayfeh-Rezaie (last updated 14 October 2026) lists quaternary Hadamard matrices
  QH590, QH638, QH782 and QH926 under "unpublished approaches". They are of this Yamada–Kiyasu kind (we recovered
  U = (I, 0, S, −S) and V = (M, N, I, I) from QH590), and nauty finds monomial matrices taking each of them to a
  symmetric quaternary matrix, so Turyn's map gives symmetric Hadamard matrices of orders 1180, 1276, 1564 and 1852 from
  them too. The page says nothing about symmetry, and the real orders 4v are not listed; we report this overlap for the
  four values v = 295, 319, 391 and 463.
- The skew inputs are classical except our skew matrices of orders 820, 836, 1036 and 1148 (2026).

The sources searched are listed in [NOVELTY.md](../NOVELTY.md).

## Reproduction and verification

```
python3 search/sym_yk.py targets           # the 34 remaining census values v = q + 2, q = 5 (mod 8), and their skew inputs
python3 search/sym_yk.py certify 295       # writes symmetric/1180/certificate.json
python3 verify/verify.py symmetric/1180    # rebuilds H and checks it
```

`verify/verify.py` (construction `symmetric-yamada-kiyasu`) rebuilds the skew-Hadamard input from its nested
certificate with its own builder (Paley's matrix `paley-skew`, Goethals–Seidel, two-block, doubling, or a certificate
of [skew](../skew)) and checks it, including Q + Qᵀ = 2I; computes Yamada's circulants over GF(q²) and checks that they
are symmetric with M² + N² = qI; builds 2X and 2Y in integers and checks that X and Y are symmetric with disjoint
supports covering every position; and checks H = Hᵀ, H Hᵀ = 4vI and the SHA-256 of the canonical CSV.

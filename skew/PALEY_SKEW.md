# Paley's matrix in sixteen blocks, times three: skew-Hadamard matrices of order 3(q + 1)

**Theorem.** For every prime power q ≡ 3 (mod 4) there is a skew-Hadamard matrix of order 3(q + 1).

Paley's skew-Hadamard matrix of order q + 1 = 4s is a 4 × 4 array of negacyclic s × s blocks. An anti-automorphism of
Paley's matrix (Frobenius) pairs two blocks of one block column, and this pairing is exactly what the T-matrices of
order 3 need to make the first Goethals–Seidel block skew-type. The block form and the pairing are Yamada's (1991; for
q ≡ 3 mod 8 it is Theorem 1 of Seberry–Yamada 1993: a skew A and B₁ = −B₀ᵀ), and the multiplication by T-matrices is
Seberry and Yamada's (1990); what we add is that the product keeps a skew-type block, which no source we found states.
The construction gives skew-Hadamard matrices of 52 orders 4v that the Đoković–Seberry census lists as undecided
(q ≡ 3 mod 8), and of six orders 8v and 16v (q ≡ 7 mod 8), five of them below the smallest power of two in
Cati–Pasechnik Table 5; the sixth, 7656 = 8 · 957, is Known (Wallis–Whiteman, see the Status). Seven more census
values in the family, 261, 285, 429, 465, 513, 663 and 2475, were already settled in this repository by searches
([index](README.md)).

| Order | v | q | s | Census; Table 5 | Status |
|------:|---|---:|---:|---|---|
| [996](996) | 249 = 3 · 83 | 331 | 83 | undecided; Table 5: 249(4) | New, obtainable classically |
| [1476](1476) | 369 = 3² · 41 | 491 | 123 | undecided; Table 5: 369(4) | New, obtainable classically |
| [1692](1692) | 423 = 3² · 47 | 563 | 141 | undecided; Table 5: 423(4) | New, obtainable classically |
| [1764](1764) | 441 = 3² · 7² | 587 | 147 | undecided; Table 5: 441(3) | New, obtainable classically |
| [1800](1800) | 225 = 3² · 5² (order 8 · v) | 599 | 150 | Table 5: 225(4) | New, obtainable classically |
| [2076](2076) | 519 = 3 · 173 | 691 | 173 | undecided; Table 5: 519(4) | New, obtainable classically |
| [2436](2436) | 609 = 3 · 7 · 29 | 811 | 203 | undecided; Table 5: 609(3) | New, obtainable classically |
| [2484](2484) | 621 = 3³ · 23 | 827 | 207 | undecided; Table 5: 621(3) | New, obtainable classically |
| [2724](2724) | 681 = 3 · 227 | 907 | 227 | undecided; Table 5: 681(4) | New, obtainable classically |
| [2916](2916) | 729 = 3⁶ | 971 | 243 | undecided; Table 5: 729(4) | New, obtainable classically |
| [3060](3060) | 765 = 3² · 5 · 17 | 1019 | 255 | undecided; Table 5: 765(4) | New, obtainable classically |
| [3156](3156) | 789 = 3 · 263 | 1051 | 263 | undecided; Table 5: 789(3) | New, obtainable classically |
| [3516](3516) | 879 = 3 · 293 | 1171 | 293 | undecided; not in Table 5 | First explicit; see its README |
| [3564](3564) | 891 = 3⁴ · 11 | 1187 | 297 | undecided; Table 5: 891(3) | New, obtainable classically |
| [3876](3876) | 969 = 3 · 17 · 19 | 1291 | 323 | undecided; Table 5: 969(4) | New, obtainable classically |
| [4344](4344) | 543 = 3 · 181 (order 8 · v) | 1447 | 362 | Table 5: 543(4) | New, obtainable classically |
| [4356](4356) | 1089 = 3² · 11² | 1451 | 363 | undecided | New, obtainable classically |
| [4380](4380) | 1095 = 3 · 5 · 73 | 1459 | 365 | undecided | New, obtainable classically |
| [4500](4500) | 1125 = 3² · 5³ | 1499 | 375 | undecided | New, obtainable classically |
| [4572](4572) | 1143 = 3² · 127 | 1523 | 381 | undecided | New, obtainable classically |
| [4596](4596) | 1149 = 3 · 383 | 1531 | 383 | undecided | New, obtainable classically |
| [4632](4632) | 579 = 3 · 193 (order 8 · v) | 1543 | 386 | Table 5: 579(5) | New, obtainable classically |
| [4740](4740) | 1185 = 3 · 5 · 79 | 1579 | 395 | undecided | New, obtainable classically |
| [4860](4860) | 1215 = 3⁵ · 5 | 1619 | 405 | undecided | New, obtainable classically |
| [4884](4884) | 1221 = 3 · 11 · 37 | 1627 | 407 | undecided | New, obtainable classically |
| [5364](5364) | 1341 = 3² · 149 | 1787 | 447 | undecided | New, obtainable classically |
| [5436](5436) | 1359 = 3² · 151 | 1811 | 453 | undecided | New, obtainable classically |
| [5604](5604) | 1401 = 3 · 467 | 1867 | 467 | undecided | New, obtainable classically |
| [5796](5796) | 1449 = 3² · 7 · 23 | 1931 | 483 | undecided | New, obtainable classically |
| [5964](5964) | 1491 = 3 · 7 · 71 | 1987 | 497 | undecided | New, obtainable classically |
| [6036](6036) | 1509 = 3 · 503 | 2011 | 503 | undecided | New, obtainable classically |
| [6084](6084) | 1521 = 3² · 13² | 2027 | 507 | undecided | New, obtainable classically |
| [6396](6396) | 1599 = 3 · 13 · 41 | 2131 | 533 | undecided | New, obtainable classically |
| [6540](6540) | 1635 = 3 · 5 · 109 | 2179 | 545 | undecided | New, obtainable classically |
| [6612](6612) | 1653 = 3 · 19 · 29 | 2203 | 551 | undecided | New, obtainable classically |
| [6732](6732) | 1683 = 3² · 11 · 17 | 2243 | 561 | undecided | New, obtainable classically |
| [6756](6756) | 1689 = 3 · 563 | 2251 | 563 | undecided | New, obtainable classically |
| [7236](7236) | 1809 = 3³ · 67 | 2411 | 603 | undecided | New, obtainable classically |
| [7272](7272) | 909 = 3² · 101 (order 8 · v) | 2423 | 606 | Table 5: 909(4) | New, obtainable classically |
| [7380](7380) | 1845 = 3² · 5 · 41 | 2459 | 615 | undecided | New, obtainable classically |
| [7620](7620) | 1905 = 3 · 5 · 127 | 2539 | 635 | undecided | New, obtainable classically |
| [7656](7656) | 957 = 3 · 11 · 29 (order 8 · v) | 2551 | 638 | Table 5: 957(4) | Known (Wallis–Whiteman) |
| [7740](7740) | 1935 = 3² · 5 · 43 | 2579 | 645 | undecided | New, obtainable classically |
| [8052](8052) | 2013 = 3 · 11 · 61 | 2683 | 671 | undecided | New, obtainable classically |
| [8100](8100) | 2025 = 3⁴ · 5² | 2699 | 675 | undecided | New, obtainable classically |
| [8196](8196) | 2049 = 3 · 683 | 2731 | 683 | undecided | New, obtainable classically |
| [8412](8412) | 2103 = 3 · 701 | 2803 | 701 | undecided | New, obtainable classically |
| [8532](8532) | 2133 = 3³ · 79 | 2843 | 711 | undecided | New, obtainable classically |
| [8556](8556) | 2139 = 3 · 23 · 31 | 2851 | 713 | undecided | New, obtainable classically |
| [8892](8892) | 2223 = 3² · 13 · 19 | 2963 | 741 | undecided | New, obtainable classically |
| [8916](8916) | 2229 = 3 · 743 | 2971 | 743 | undecided | New, obtainable classically |
| [9252](9252) | 2313 = 3² · 257 | 3083 | 771 | undecided | New, obtainable classically |
| [9564](9564) | 2391 = 3 · 797 | 3187 | 797 | undecided | New, obtainable classically |
| [9612](9612) | 2403 = 3³ · 89 | 3203 | 801 | undecided | New, obtainable classically |
| [9780](9780) | 2445 = 3 · 5 · 163 | 3259 | 815 | undecided | New, obtainable classically |
| [9972](9972) | 2493 = 3² · 277 | 3323 | 831 | undecided | New, obtainable classically |
| [9996](9996) | 2499 = 3 · 7² · 17 | 3331 | 833 | undecided | New, obtainable classically |
| [14352](14352) | 897 = 3 · 13 · 23 (order 16 · v) | 4783 | 1196 | Table 5: 897(5) | New, obtainable classically |

## Construction

**1. Paley's matrix in closed form.** Let q ≡ 3 (mod 4) be a prime power, so −1 is a non-square in GF(q) and
GF(q²) = GF(q)(i) with i² = −1. Write ū = u^q and Im(u₁ + u₂i) = u₂. Let z be a primitive element of GF(q²) and
a = z^(q+1), a generator of GF(q)*, so χ(a) = −1 (χ the quadratic character of GF(q)). The q + 1 points of PG(1, q)
are the classes of x_m = z^m, m ∈ ℤ_{q+1}, and the symplectic form is det(u, v) = Im(ū v). Since
det(x_m, x_n) = Im(z^{qm+n}) = a^m Im(z^{n−m}), Paley's skew-Hadamard matrix, on any choice of representatives
x_m = z^m (m running over a set of residues mod q + 1), is

    P(x_m, x_n) = [m = n] + (−1)^m F(n − m),      F(d) = χ(Im z^d),

and different choices of representatives change P only by a diagonal similarity D P D (D = diag(±1)). F has two
symmetries:

- F(d + q + 1) = χ(a) F(d) = −F(d), since z^(q+1) = a ∈ GF(q);
- F(−d) = (−1)^(d+1) F(d), since z^(−d) = conj(z^d) / a^d and Im ū = −Im u (this is Frobenius, an anti-automorphism
  of P: det(ū, v̄) = −det(u, v)).

**2. Sixteen negacyclic blocks.** Let s = (q + 1)/4. Multiplication by z⁴ has the four orbits O_k = {x_m : m ≡ k
(mod 4)} of length s. List O_k as x_{k'+4i}, i = 0, …, s − 1, with base exponents k' = 0, 1, 2 and −1 (for O₃, the
base point x_{−1} = x_q up to the sign χ(a)). Then the block of P in rows O_k and columns O_j is negacyclic:

    N_kj[i, i'] = Ñ_kj(i' − i),    Ñ_kj(y) = [k = j, y = 0] + (−1)^k F(j' − k' + 4y),    Ñ_kj(y + s) = −Ñ_kj(y).

Because P + Pᵀ = 2I, the diagonal blocks are skew-type (N_kk + N_kkᵀ = 2I), and the block rows are orthogonal, Σ_j N_kj
N_ljᵀ = 4s δ_kl I. (For s odd, negacyclic and circulant blocks are equivalent;
[PALEY_BLOCKS.md](../hadamard/PALEY_BLOCKS.md) builds circulant blocks through a torus element of PSL(2, q). For s even,
every lift of that element to SL(2, q) has order 2s and the blocks are necessarily negacyclic.)

**3. The pairing.** In block column 0, by the second symmetry of F with d = 1 + 4y,

    Ñ₃₀(y) = −F(1 + 4y) = −F(−1 − 4y) = Ñ₁₀(−y),    so N₃₀ = N₁₀ᵀ.

(Frobenius fixes the orbits O₀ and O₂ and swaps O₁ and O₃; the base exponent −1 of O₃ absorbs the shift.)

**4. Times the T-matrices of order 3.** Let Z be the cyclic shift of order 3 and put T₀ = I, T₁ = Z, T₂ = 0,
T₃ = −Z²: circulant (0, ±1) matrices with supports {0}, {1}, ∅, {2} partitioning ℤ₃ and Σ T_k T_kᵀ = 3I. Then

    A_j = Σ_k T_k ⊗ N_kj        (j = 0, …, 3)

are ±1 matrices in the commutative algebra of circulant ⊗ negacyclic matrices on ℤ₃ × ℤ_s, and
Σ_j A_j A_jᵀ = Σ_{k,l} T_k T_lᵀ ⊗ Σ_j N_kj N_ljᵀ = 12s I (Seberry–Yamada's multiplication theorem; see
[PALEY_BLOCKS.md](../hadamard/PALEY_BLOCKS.md)). The first one is skew-type:

    A₀ + A₀ᵀ = I ⊗ (N₀₀ + N₀₀ᵀ) + (Z ⊗ N₁₀ + Z² ⊗ N₁₀ᵀ) − (Z² ⊗ N₁₀ᵀ + Z ⊗ N₁₀) = 2I.

**5. The Goethals–Seidel array.** With R = R₃ ⊗ J_s (R₃ : x ↦ −x on ℤ₃, J_s the back-identity), R X R = Xᵀ for every X
in the algebra (J N J = Nᵀ holds for negacyclic N as for circulant N), so

```
    [  A₀    A₁R    A₂R    A₃R  ]
H = [ -A₁R   A₀     A₃ᵀR  -A₂ᵀR ]
    [ -A₂R  -A₃ᵀR   A₀     A₁ᵀR ]
    [ -A₃R   A₂ᵀR  -A₁ᵀR   A₀   ]
```

is a Hadamard matrix of order 12s = 3(q + 1), and H + Hᵀ = 2I because A₀ + A₀ᵀ = 2I and (A_jR)ᵀ = A_jR. ∎

`search/skew_3q.py q` builds H from these formulas for a prime q and checks it; it passes for all 86 primes
q ≡ 3 (mod 4) from 7 to 991. The certificates use the coordinates of [PALEY_BLOCKS.md](../hadamard/PALEY_BLOCKS.md)
instead (Paley's matrix on PG(1, q), a torus element of order s, prime powers included). There `search/paley_blocks.py
skew q` finds by search a block column c, rows b, d, a shift δ and a sign ε with N_dc(y) = ε N_bc(δ − y) (always
c = 0), moves the base point of O_d, and puts T_c = I, T_b = Z, T_d = −εZ², 0 on the fourth orbit. For q ≡ 7 (mod 8)
its blocks are negacyclic and the certificate says so.

**What does not extend.**

- Other multipliers. In general the first block A_c = Σ_k T_k ⊗ N_kc is skew-type only if T_c = I and the other three
  are X, −εXᵀ and a skew Y (on the orbit whose block is symmetric up to a shift), with 2 PAF(X) + PAF(Y) = (t − 1)δ.
  The trivial character forces t = 2a² + 1 (Y has sum 0). An exhaustive search over ℤ₉, ℤ₃ × ℤ₃ and ℤ₁₉ found no such
  structure; t = 3 (X = Z, Y = 0) is the only one we have.
- A symmetric analogue. The same blocks give a symmetric A_c (T_c = 0, T_o = I on the other fixed orbit), but the
  propus array also needs two of the other three blocks with equal Gram matrices (A_j A_jᵀ = A_k A_kᵀ; one of them is
  then used twice). For q = 19, 43, 67 no two of them are related by any affine map of ℤ₃ × ℤ_s, and an exhaustive
  check of all 192 placements of the T-matrices of order 3 and all base points of the orbits finds equal Gram matrices
  for q = 11 only (order 36), none for q = 7, 19, 23, 31 and 43.

## Status

**New, but obtainable classically** for 51 orders 4v and five orders 8v and 16v; **First explicit construction** for
3516 = 4 · 879; **Known** for 7656 = 8 · 957. Ordinary Hadamard matrices of all these orders are obtainable classically
(the multiplication without the pairing; for 3156, 6036 and 6756 this repository built the first ones,
[PALEY_BLOCKS.md](../hadamard/PALEY_BLOCKS.md)); what is new is the skew property. Every ingredient is published: the
sixteen-block form with the transposed pair (Yamada 1991, restated as Theorem 1 of Seberry–Yamada, JCMCC 13, 1993,
whose Theorem 6 multiplies such matrices by T-matrices under amicability conditions that T-matrices of order 3 fail),
the multiplication (Seberry–Yamada 1990, Theorem 4) and the Goethals–Seidel criterion for skew matrices. Balonin and
Đoković (2015, §9) noticed for q = 139 that one block of the array is skew-type and that two others are reverses of
each other. We found no source that combines these into a skew matrix of order 3(q + 1).

- **7656.** Wallis and Whiteman (Bull. Austral. Math. Soc. 7, 1972) give skew-Hadamard matrices of order 4(q + 1) for
  prime powers q ≡ 9 (mod 16) (Momihara–Xiang, arXiv:1801.08776, Propositions 1.3–1.4 and Corollary 1.7), and 7656 =
  4(1913 + 1). Cati–Pasechnik's 957(4) misses this; the matrix is built and checked here. 14352 = 16 · 897 also follows
  from Đoković's skew matrix of order 276 times the amicable pair of order 52 of [2652](2652). 3060: Table 5's 765(4)
  misses the family (|t| + 1)(q + 1) with q = 509, which gives 8 · 765; 4 · 765 is undecided in the census.

- **The census and the tables.** The Đoković–Seberry census of skew orders 4v (odd v < 2500, May 2025) lists only the
  Paley, Szekeres, Szekeres–Whiteman and Spence series and 41 sporadic values, and all 59 values v = 3(q + 1)/4 of the
  family (q ≡ 3 mod 8) as undecided. Cati–Pasechnik Table 5 gives t ≥ 3 for the 18 of them below 1000 except 663 and
  879, and t = 4 or 5 for the six values of the q ≡ 7 (mod 8) orders. Seberry–Yamada (1992), Table 7.1, Seberry (1978),
  Seberry–Yamada (1991) and Table 1.51 of the Handbook (2007) agree for the values they cover. 663 and 879 are blank
  (t = 2) in Seberry–Yamada's Table 7.1, and hence absent from Cati–Pasechnik's Table 5, with no construction cited; we
  read these blanks as errors (also 783; for all three, 2v − 1 ≡ 5 (mod 8) is composite, which suggests a misapplied
  Szekeres entry), so 3516 is labelled First explicit construction, as [2652](2652) is.
- **The families.** None of the skew families listed by Seberry (1978), Seberry–Yamada (1992, SI–SXI), the Handbook,
  Koukouvinos–Stylianou (2008, first pages) or the census gives order 3(q + 1) except when an unrelated number happens
  to qualify (3q + 2 a prime power for Paley, (3q + 1)/2 for Szekeres), which happens for none of the 58 orders here.
  Seberry's Theorem 7 (1978) with n = 4 would give 3(q + 1) from an orthogonal design OD(q + 1; 1, (q + 1)/4,
  3(q + 1)/4 − 1); we found no source that constructs these designs.
- **The ingredients.** The four-orbit structure of Paley's matrix is classical (Whiteman 1972 for the Goethals–Seidel
  form, as summarized by Seberry–Yamada 1992; Yamada 1991 and Seberry–Yamada 1993 for the generalized quaternion type,
  with the transposed pair;
  Ito's dicyclic Hadamard groups; Xia–Xia–Seberry 2005 for an infinite family of Goethals–Seidel arrays at order q + 1),
  and so is the multiplication by T-matrices (Seberry–Yamada 1990, Theorem 4; 1992, Theorem 3.10), which gives an
  ordinary Hadamard matrix of order 3(q + 1). Xia, Xia, Seberry and Wu (2006) multiply negacyclic Goethals–Seidel
  matrices by base sequences, also for ordinary Hadamard matrices. Xia–Xia–Seberry (2005, Theorem 3) use the
  conjugation x ↦ x^q, the anti-automorphism behind the pairing. None of these sources states a skew conclusion.
- **Recent work and software.** Đoković (arXiv:2508.17141, 2025) uses the Xia–Xia–Seberry blocks only at order
  q + 1. The arXiv skew-Hadamard records, the Hadamard records of 2025–2026, Zenodo, the JD-Jones-ASES repositories and
  Schumer's deposit have none of the 58 orders. SageMath's development version lists 996 in `unknown_skew_hadamard`,
  and its skew constructions give none of the 58 orders.
- **Not read.** Whiteman (Discrete Math. 2, 1972); Seberry–Yamada's book (2020), chapter 9 and its Table 9.2;
  Haviland's thesis (BYU, 2023) and Anderson–Haviland–Holmes–Humphries–Magland (Graphs Combin. 39, 2023) on skew
  relative difference sets in dicyclic groups, whose abstract describes one infinite family built from Ito's and
  Schmidt's results. Our matrices are not group-developed in that way: for q = 19, 43 and 59 their automorphism
  groups have order 2. The sources searched are listed in [NOVELTY.md](../NOVELTY.md).

## Reproduction and verification

    python3 search/paley_blocks.py skew 331            # writes skew/996/certificate.json and sequences.txt
    python3 verify/verify.py skew/996                  # rebuilds and checks H Hᵀ = 996 I and H + Hᵀ = 2I
    python3 search/skew_3q.py 331 599                  # the closed formula, for primes q ≡ 3 (mod 4)

`verify/verify.py` (construction `paley-blocks` with `"skew": true`) builds Paley's matrix independently from the
symplectic form on PG(1, q), finds the blocks from the certificate's τ, power and c (negacyclic when s is even), moves
the base points by the certificate's shifts, takes the block columns in its order, rebuilds the four sequences from the
T-matrices and compares them with the certificate. It then checks the sum of the (twisted) periodic autocorrelations,
assembles the Goethals–Seidel matrix, and checks H Hᵀ = nI, H + Hᵀ = 2I and the SHA-256 of the canonical CSV. All 58
certificates pass. The matrices are not stored.

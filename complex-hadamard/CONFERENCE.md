# Quaternary conference matrices and the Yamada–Kiyasu template

The quaternary Hadamard matrices of orders 1286, 1718, 2102 and 3398, and by Turyn's doubling the Hadamard matrices of
orders 2572, 3436, 4204 and 6796, come from quaternary conference matrices of orders 322, 430, 526 and 850 by one
template. This file states and proves its three parts: an array that turns four sequences on an abelian group into a
quaternary conference matrix (plain and bordered), the template that turns a quaternary conference matrix of order
(q + 3)/2 into a quaternary Hadamard matrix of order 2(q + 2), and the doubling. `verify/verify.py` rebuilds every
step from the certificate (construction `complex-yamada-kiyasu`).

| order 2N | N | q | core order (q + 3)/2 | core | symmetry group of the sequences | Hadamard matrix |
|---:|---|---:|---:|---|---|---|
| 1286 | 643 | 641 | 322 | zero-diagonal array over ℤ₁₆₁ = F₇ × ℤ₂₃ | ⟨93⟩, order 3 | [2572](../hadamard/2572) |
| 1718 | 859 | 857 | 430 | zero-diagonal array over ℤ₂₁₅ | ⟨49⟩, order 6 | [3436](../hadamard/3436) |
| 2102 | 1051 | 1049 | 526 | bordered zero-diagonal array over ℤ₂₆₂ | ⟨53⟩, order 5 | [4204](../hadamard/4204) |
| 3398 | 1699 | 1697 | 850 | zero-diagonal array over ℤ₅² × ℤ₁₇ | order 24 | [6796](../hadamard/6796) |

Labels: 1286 and 1718 are New; 2102 and 3398 are New, beyond the tables ([NOVELTY.md](../NOVELTY.md)). Yamada's
theorem needs a symmetric conference matrix, which cannot exist for these cores (321, 429, 525 and 849 are not sums of
two squares); the template needs no symmetry.

## Notation

- A **unit** is one of 1, −1, i, −i.
- A **quaternary conference matrix** of order c is a c × c matrix Q with zero diagonal, entries 1, −1, i, −i off the
  diagonal and Q Q* = (c − 1) I. Then also Q* Q = (c − 1) I, since Q is a multiple of a unitary matrix. Real symmetric
  and skew conference matrices are examples. Scaling rows and columns by units keeps these properties, so Q can be
  brought to the **normal form** Q = [[0, e], [eᵀ, Z]] with e the all-one row of length n = c − 1. Then Z has zero
  diagonal and entries 1, −1, i, −i elsewhere, and comparing Q Q* and Q* Q with (c − 1) I gives Z Z* = Z* Z = nI − J
  and zero row and column sums of Z (J is the all-one matrix).
- G is an abelian group of order m. For f : G → ℂ let X_f[g, h] = f(h − g) and PAF_f(s) = Σ_t f(t) f(t + s). Group
  matrices commute, (X_f X_fᵀ)[g, h] = PAF_f(g − h), and the row and column sums of X_f are r_f = Σ_g f(g).
- R[g, h] = 1 if g + h = c₀ and 0 otherwise, for a fixed c₀ ∈ G (c₀ = 0 unless the certificate says otherwise). Then
  R = Rᵀ, R² = I and R X_f R = X_fᵀ for every f, hence R X_fᵀ = X_f R.
- ω = (1 + i)/2, so ω ω̄ = 1/2, ω² = i/2, ω̄² = −i/2. For x, y = ±1, ωx + ω̄y is one of 1, i, −i, −1.
- P = [[1, −1], [−1, 1]] and J₂ is the all-one matrix of order 2: P J₂ = J₂ P = 0, P² = 2P, J₂² = 2J₂, P + J₂ = 2I.

## 1. The zero-diagonal array

**Theorem 1.** Let A, B, C, D : G → {0, ±1} with A, B symmetric (A(−g) = A(g), B(−g) = B(g)), A(0) = B(0) = 0, A and B
equal to ±1 elsewhere, C and D equal to ±1 everywhere, and

    PAF_A(g) + PAF_B(g) + PAF_C(g) + PAF_D(g) = 0   for g ≠ 0.

(At g = 0 the sum is 4m − 2.) Put X = ω X_A + ω̄ X_B and Y = ω X_C + ω̄ X_D R. Then

    Q = [ X    Y  ]
        [ Y*  −X* ]

is a quaternary conference matrix of order 2m.

*Proof.* The entry X[g, h] is 0 for g = h and one of ±1, ±i otherwise; every entry of Y is a unit. So Q has zero
diagonal and entries ±1, ±i elsewhere.

X is symmetric, and X_A, X_B commute with R because R X_A R = X_Aᵀ = X_A. Hence X Y = Y X.

Since X_A, X_B are real, symmetric and commute,

    X X* = ½(X_A X_Aᵀ + X_B X_Bᵀ) + (i/2)(X_A X_Bᵀ − X_B X_Aᵀ) = ½(X_A X_Aᵀ + X_B X_Bᵀ).

With D′ = X_D R we have D′D′ᵀ = X_D X_Dᵀ and X_C D′ᵀ = X_C R X_Dᵀ = X_C X_D R = X_D X_C R = D′ X_Cᵀ, so in the same
way

    Y Y* = ½(X_C X_Cᵀ + X_D X_Dᵀ).

Therefore X X* + Y Y* = ½ Σ_f X_f X_fᵀ, whose (g, h) entry is ½ Σ_f PAF_f(g − h). This is ½(4m − 2) = 2m − 1 for g = h
and 0 otherwise. X\*X and Y\*Y expand to the same expressions, because all group matrices (transposes included)
commute and R X_fᵀ R = X_f, so X\*X + Y\*Y = (2m − 1) I as well. Then

    Q Q* = [ X X* + Y Y*     X Y − Y X     ] = (2m − 1) I.   ∎
           [ (X Y − Y X)*    X*X + Y*Y     ]

The certificates of 3436 and 6796 hold such quadruples over ℤ₂₁₅ and over ℤ₅² × ℤ₁₇.

## 2. The bordered zero-diagonal array

**Theorem 2.** Let G have even order m and let A, B, C, D be as in Theorem 1, but with

    PAF_A(g) + PAF_B(g) + PAF_C(g) + PAF_D(g) = −4   for g ≠ 0.

Then the row sums are r_A, r_B = ±1 and r_C = r_D = 0. Put ρ = ω r_A + ω̄ r_B, a unit. Let t and v be units
(elements of {1, −1, i, −i}) and

    s = ρ² t̄,   y = −ρ t̄ v,   u = ρ² t̄² v,
    X′ = ω X_A + ω̄ X_B,   Y′ = ω X_C + ω̄ X_D R,

    X = [ 0     s e ]      Y = [ y     u e ]      Q = [ X    Y  ]
        [ t eᵀ  X′  ]          [ v eᵀ  Y′  ]          [ Y*  −X* ]

with e the all-one row of length m. Then Q is a quaternary conference matrix of order 2(m + 1).

*Proof.* Summing the condition over g gives Σ_f r_f² = (4m − 2) − 4(m − 1) = 2. A and B have m − 1 entries ±1, an odd
number because m is even, so r_A and r_B are odd; hence r_A² = r_B² = 1 and r_C = r_D = 0.

As in Theorem 1, X′ Y′ = Y′ X′, X′ and Y′ are normal, and X′X′* + Y′Y′* = ½ Σ_f X_f X_fᵀ = (2m + 1) I − 2J (the
entries of Σ_f X_f X_fᵀ are 4m − 2 on the diagonal and −4 elsewhere). Since X_f eᵀ = r_f eᵀ and e X_f = r_f e (and
e R = e),

    X′ eᵀ = ρ eᵀ,   e X′ = ρ e,   Y′ eᵀ = 0,   e Y′ = 0.

Multiplying out the blocks, with e eᵀ = m and eᵀe = J:

    X X* + Y Y* = [ m + (1 + m)         (s ρ̄ + y v̄) e   ]
                  [ (s̄ ρ + v ȳ) eᵀ     2J + X′X′* + Y′Y′* ],
    X*X + Y*Y   = [ m + (1 + m)         (t̄ ρ + ȳ u) e   ]
                  [ (t ρ̄ + ū y) eᵀ     2J + X′*X′ + Y′*Y′ ],
    X Y − Y X   = [ (s v − u t) m       −(y s + u ρ) e   ]
                  [ (t y + v ρ) eᵀ      (t u − v s) J    ].

The lower right blocks of the first two are (2m + 1) I. The remaining conditions

    s ρ̄ + y v̄ = 0,   t̄ ρ + ȳ u = 0,   y s + u ρ = 0,   t y + v ρ = 0,   s v = u t

hold for the stated s, y, u, because |ρ| = |t| = |v| = 1: for example s ρ̄ + y v̄ = ρ t̄ − ρ t̄ = 0,
ȳ u = −ρ t̄, y s = −ρ³ t̄² v = −u ρ, t y = −ρ v and s v = ρ² t̄ v = u t. So Q Q* = (2m + 1) I. Zero diagonal: the
corner of X is 0, the diagonal of X′ is 0, and every other entry of Q is a unit. ∎

The group, border and sequences of 4204 are those of the case m = 262, r_A = r_B = 1, ρ = 1, t = v = 1, which gives
s = 1, y = −1, u = 1. Theorem 2 relates to Theorem 1 as the bordered Szöllősi array of [BORDERED.md](BORDERED.md)
relates to Szöllősi's array. Theorem 1 is Szöllősi's array ([README.md](README.md#construction)) over an abelian group
with A(0) = B(0) = 0, so that the output has zero diagonal; Theorem 2 adds a border.

## 3. The template

**Theorem 3.** Let q be a prime, q ≡ 1 (mod 4), n = (q + 1)/2 (odd) and v = q + 2 = 2n + 1. Let Q = [[0, e], [eᵀ, Z]]
be a quaternary conference matrix of order n + 1 = (q + 3)/2 in normal form. Let M and N be symmetric circulant
matrices of order n with M of zero diagonal and ±1 elsewhere, N with all entries ±1, and M² + N² = qI (Yamada's
circulants, below). With I the identity of order n and e₂ the all-one row of length 2n, put

    A = P ⊗ M + J₂ ⊗ I,   B = P ⊗ N,   X′ = ω A + ω̄ B,   Y′ = J₂ ⊗ Z + P ⊗ I,

    X = [ 1     e₂ ]      Y = [ 1      i e₂ ]      K = [ X    Y  ]
        [ −e₂ᵀ  X′ ]          [ −i e₂ᵀ  Y′  ]          [ Y*  −X* ].

Then K has entries 1, −1, i, −i and K K* = 2v I: a quaternary Hadamard matrix of order 2(q + 2). No symmetry of Q is
used.

**Yamada's circulants.** Let ψ be the quadratic character of GF(q) with ψ(0) = 0, r a non-square modulo q and ξ a
primitive element of GF(q²) = GF(q)[w]/(w² − r), with trace Tr(a + bw) = 2a. M and N are the circulants with first
rows

    β_m = ψ(Tr ξ^(4m + n))  (M),   α_m = ψ(Tr ξ^(4m))  (N),   0 ≤ m < n.

They are symmetric and satisfy M N = N M and M² + N² = qI (Yamada's Lemma 3; M. Yamada, Some new series of Hadamard
matrices, J. Austral. Math. Soc. Ser. A 46 (1989) 371–383, as used in
[symmetric/YAMADA_KIYASU.md](../symmetric/YAMADA_KIYASU.md) and [hadamard/15292](../hadamard/15292)). For each
certificate the builder checks the symmetry, the entry pattern and PAF_α + PAF_β = q δ, which is M² + N² = qI. The
pattern is elementary: Tr ξ^k = 0 exactly when ξ^(k(q−1)) = −1, that is k ≡ n (mod 2n); as n is odd, 4m + n ≡ n (mod
2n) only for m = 0, and 4m ≡ n (mod 2n) never. So β_0 = 0, all other β_m and all α_m are ±1.

*Proof.* **Entries.** The entries of A and B at ((j, a), (k, b)) are P_ab M_jk + δ_jk and P_ab N_jk. For j ≠ k they
are ±1, and for j = k they are 1 and ±N_jj = ±1 (M_jj = 0). So every entry of X′ is a unit. Y′ has the entry Z_jk, a
unit, for j ≠ k and P_ab = ±1 for j = k (Z_jj = 0). The border entries are units.

**Row sums.** (P ⊗ M)(1) = (P1) ⊗ (M1) = 0 and (J₂ ⊗ I)(1) = 2 · 1, so A eᵀ = 2 eᵀ and B eᵀ = 0 (1 is the all-one
column). A and B are symmetric, so e₂ A = 2 e₂ and e₂ B = 0. Z has zero row and column sums, so Y′ e₂ᵀ = 0 and
e₂ Y′ = 0. Hence X′ e₂ᵀ = (1 + i) e₂ᵀ and e₂ X′ = (1 + i) e₂.

**The cores.** Using P J₂ = J₂ P = 0, P² = 2P and J₂² = 2J₂:

    A Y′ = 2 P ⊗ M + 2 J₂ ⊗ Z = Y′ A,   B Y′ = 2 P ⊗ N = Y′ B,   so X′ Y′ = Y′ X′;
    A B = B A = 2 P ⊗ M N,   A² = 2 P ⊗ M² + 2 J₂ ⊗ I,   B² = 2 P ⊗ N².

A and B are real, symmetric and commute, so X′ is normal and
X′X′* = ½(A² + B²) = P ⊗ (M² + N²) + J₂ ⊗ I = q P ⊗ I + J₂ ⊗ I. Next Y′Y′* = Y′\*Y′ = 2 J₂ ⊗ Z Z* + 2 P ⊗ I = 2 J₂ ⊗
(nI − J_n) + 2 P ⊗ I. With q + 2 = v, J = J₂ ⊗ J_n and P + J₂ = 2I,

    X′X′* + Y′Y′* = v (P + J₂) ⊗ I − 2 J₂ ⊗ J_n = 2v I − 2J,

and the same for X′\*X′ + Y′\*Y′.

**The borders.** With e₂ e₂ᵀ = 2n, e₂ᵀe₂ = J and the row sums above, the blocks of X X* + Y Y* are

    (1, 1): (1 + 2n) + (1 + 2n) = 2v;
    (1, 2): −e₂ + e₂ X′* + i e₂ + i e₂ Y′* = −e₂ + (1 − i) e₂ + i e₂ + 0 = 0;
    (2, 2): J + X′X′* + J + Y′Y′* = 2J + 2v I − 2J = 2v I,

and (2, 1) is the adjoint of (1, 2). For X\*X + Y\*Y the (1, 1) block is again 2v, the (1, 2) block is
(e₂ − e₂ X′) + (i e₂ + i e₂ Y′) = e₂ − (1 + i) e₂ + i e₂ = 0, and the (2, 2) block is 2J + 2v I − 2J. For X Y and Y X
the (1, 1) blocks are both 1 − 2ni, the (1, 2) blocks are both i e₂ (e₂ Y′ = 0 and e₂ + i e₂ X′ = i e₂), the (2, 1)
blocks are both −i e₂ᵀ, and the (2, 2) blocks are −iJ + X′Y′ and −iJ + Y′X′, which agree. So
X X* + Y Y* = X\*X + Y\*Y = 2v I and X Y = Y X, and K K* = 2v I. ∎

The certificates build Q from the sequences of Theorem 1 or 2 (permuting columns to put the zeros on the diagonal if
necessary, then scaling to normal form), check Q Q* = (c − 1) I and M² + N² = qI, and assemble K as above.

## 4. Turyn's doubling

**Theorem 4.** If K = K_r + i K_i is a quaternary Hadamard matrix of order 2v, then

    H = [  K_r + K_i    K_r − K_i ]
        [ −(K_r − K_i)  K_r + K_i ]

is a Hadamard matrix of order 4v.

*Proof.* In every position exactly one of K_r, K_i is nonzero, so X = K_r + K_i and Y = K_r − K_i have entries ±1.
K K* = 2v I gives K_r K_rᵀ + K_i K_iᵀ = 2v I and K_i K_rᵀ = K_r K_iᵀ. Then X Xᵀ + Y Yᵀ = 2(K_r K_rᵀ + K_i K_iᵀ) = 4v I
and Y Xᵀ − X Yᵀ = 2(K_r K_iᵀ − K_i K_rᵀ) = 0, so H Hᵀ = [[X Xᵀ + Y Yᵀ, Y Xᵀ − X Yᵀ], [X Yᵀ − Y Xᵀ, Y Yᵀ + X Xᵀ]] = 4v
I. ∎

With v = q + 2 this gives the orders 4(q + 2) = 2572, 3436, 4204 and 6796.

## Certificates and verification

The certificates of 1286, 1718, 2102 and 3398 ([C1286.json](C1286.json), [C1718.json](C1718.json),
[C2102.json](C2102.json), [C3398.json](C3398.json)) and of 2572, 3436, 4204 and 6796 (`certificate.json` in
[hadamard](../hadamard)) have construction `complex-yamada-kiyasu` and these fields:

- `order`: 2(q + 2) in `C<n>.json` (which also has `N` = order/2), 4(q + 2) in the real certificates;
- `q`, `r`, `xi`: the prime q, the non-square r, and ξ = a + bw as `[a, b]`;
- `core`: `array` (a label, starting with "bordered" exactly when `border` is present); `group` (a string such as
  `"215"` for ℤ₂₁₅ or `"5^2x17"` for ℤ₅² × ℤ₁₇, elements indexed in mixed radix with the last factor fastest);
  `sequences` (A, B, C, D as strings over `+`, `-`, `0`, indexed by the group); optionally `reflection` (the element
  c₀ of R, as its index), `border` (`{"t": [a, b], "v": [c, d]}` for the units t = a + bi, v = c + di), `rowsums`,
  `symmetry_generators` (integer matrices acting on the coordinates of G, one per generator) and `symmetry_order`; the
  sequences must be constant on the orbits of the generators, and the row sums and the order of the group they
  generate are checked against the stated values;
- `sha256`: of the canonical text of K (complex certificates) or of H (real ones).

```
python3 verify/verify_complex.py 1286 1718 2102 3398
python3 verify/verify.py 2572 3436 4204 6796
```

`verify_complex.py` computes K K* from the real and imaginary parts with float64 matrix products, which are exact:
every entry is 0 or ±1 and every partial sum is an integer of absolute value at most 2N.

## Status and credits

The array [[X, Y], [Y*, −X*]] with A, B symmetric circulants is Szöllősi's (arXiv:2603.09572, Theorem 4, a
modification of Kharaghani and Seberry's block construction). Its zero-diagonal form over an abelian group, which
outputs conference matrices, and its bordered form (Theorems 1 and 2) are ours ([search/NOTES.md](../search/NOTES.md)
§74–75). The template is the form of the quaternary matrices QH470, QH566, QH806, QH1142 and QH1526 on [V. Zaitsev's
page](https://www.cs.uleth.ca/~zaitsevv/quaternary.html) (with Kharaghani, Mohammadian and Tayfeh-Rezaie; "Zaitsev et
al."), which states no method; we found it stated nowhere else. With a skew-Hadamard core (q ≡ 5 (mod 8)) it is built
from the data of Yamada's theorem ([symmetric/YAMADA_KIYASU.md](../symmetric/YAMADA_KIYASU.md)), and Zaitsev et al.'s
QH590, QH638, QH782, QH926, QH1310, QH1358, QH1646 and QH1886 are of that kind ([NOTES.md](../NOTES.md) §31).
Miyamoto's Corollary 4 (1991, as transcribed in the JD-Jones-ASES Hadamard-M note) reaches the same orders from a
symmetric quaternary conference matrix. The cores here have Qᵀ ≠ Q as built.

The template reaches further orders once a core is found: the FRONTIER orders 4n for n = 955, 1291, 1879, 2083, 2215,
2383, 2659, 2731, 2803, 2899 and 2971 (cores of order (n + 1)/2 = 478 to 1486; for n = 1879, 2215 and 2383 (q ≡ 5
(mod 8)) a skew-Hadamard core of order 940, 1108 or 1192 would also serve; n = 955 also gives the open complex
order 1910) and the open complex orders 1046, 1190 and 1958 (cores of order 262, 298 and 490). None of these cores
has been found ([NOTES.md](../NOTES.md) §34, [search/NOTES.md](../search/NOTES.md) §75–83).

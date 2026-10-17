# Symmetric Hadamard matrices of order 4q from any Hadamard matrix of order (q − 1)/2

**Theorem.** Let q ≡ 1 (mod 8) be a prime power. If a Hadamard matrix of order (q − 1)/2 exists, then a symmetric
Hadamard matrix of order 4q exists.

[PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md) proves this with a *symmetric* Hadamard matrix of order (q − 1)/2, from
Seberry and Yamada's Corollary 29(ii) in Seberry and Balonin's propus array. Here the input need not be symmetric. The
conditions of Seberry and Yamada's Lemma 26 split into a U-part and a V-part, and in the propus array the V-part enters
only through its square and its support. So the V-part of Corollary 29(ii), the lift of the propus array of
(I, 0, M, M), may be replaced by any symmetric weighing matrix of order 2(q − 1) and weight q with the same support, and
such a weighing matrix is built from any Hadamard matrix of order (q − 1)/2. Ordinary Hadamard matrices of order 4q
under the same hypothesis are Yamada's Theorem 3 (J. Austral. Math. Soc. Ser. A 46 (1989) 371–383).

The theorem settles the ten values of the Đoković–Seberry census (2025-05-28) that [PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md#reach)
left open, the prime powers q ≡ 1 (mod 8) for which no symmetric Hadamard matrix of order (q − 1)/2 = 4m is known
(m = 101, 107, 119, 149, 179, 239, 249, 269, 299, 305). With the 48 orders there, every census value q ≡ 1 (mod 8) that
is a prime power is now decided.

| Order | q | (q − 1)/2 | K (a Hadamard matrix of order (q − 1)/2) | Status |
|------:|---:|---:|---|---|
| [3236](3236) | 809 | 404 | Cooper–Wallis: T-matrices of order 101, Williamson matrices of order 1 | New |
| [3428](3428) | 857 | 428 | Cooper–Wallis: T-matrices of order 107, Williamson matrices of order 1 | New |
| [3812](3812) | 953 | 476 | Cooper–Wallis: T-matrices of order 17, Turyn's Williamson matrices of order 7 | New |
| [4772](4772) | 1193 | 596 | Miyamoto (q′ = 149) with T-matrices of order 37 (order 148) | New |
| [5732](5732) | 1433 | 716 | Alpöge–Voinov–Reynolds-Haertle (2026), a bordered Goethals–Seidel array over ℤ₁₇₈ | New |
| [7652](7652) | 1913 | 956 | Đoković–Kotsireas's skew Goethals–Seidel matrix (arXiv:1511.08734) | New |
| [7972](7972) | 1993 | 996 | our skew-Hadamard matrix of order 996 ([skew/996](../skew/996)) | New |
| [8612](8612) | 2153 | 1076 | Miyamoto (q′ = 269) with T-matrices of order 67 (order 268) | New |
| [9572](9572) | 2393 | 1196 | Cooper–Wallis: T-matrices of order 23, Turyn's Williamson matrices of order 13 | New |
| [9764](9764) | 2441 | 1220 | Cooper–Wallis: T-matrices of order 5, Turyn's Williamson matrices of order 61 | New |

## Proof

Let n = (q − 1)/2 (so 4 | n), m = 2n + 1 = q, J₂ the all-one matrix of order 2, D₂ = [[1, −1], [−1, 1]] and e the
all-one row of length 2n. With g a primitive element of GF(q), χ the quadratic character, s_i = g²ⁱ and t_j = g¹⁻²ʲ
(i, j < n), let

```
U₁ = −[χ(s_i − s_j)],   U₂ = [χ(s_i − t_j)],
W₁ = [[1, −e], [−eᵀ, U₁ ⊗ J₂]],   W₂ = [[1, e], [eᵀ, U₂ ⊗ J₂]],   W₃ = [[1, e], [eᵀ, 0]],
H₀ = [[W₁, W₃, W₃, W₂], [W₃, W₂, −W₁, −W₃], [W₃, −W₁, −W₂, W₃], [W₂, −W₃, W₃, −W₁]].
```

H₀ is the matrix of [PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md) with its V-part removed. Index the rows of block k
(k = 0..3) by the border k·m and the core coordinates (k, i, a) = k·m + 1 + 2i + a (i < n, a ∈ {0, 1}).

**1. The U-part.** U₁ is a symmetric circulant with zero diagonal and no other zero, U₂ is back-circulant (hence
symmetric) with no zero, U₁U₂ = U₂U₁, U₁ has row sum 1, U₂ row sum 0, and U₁² + U₂² = qI − 2J
([PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md#proof), conditions (i), (iv), (v); none involves the input). Then W₁, W₂, W₃
are symmetric and pairwise commuting (W₁W₃ = W₃W₁ needs e(U₁ ⊗ J₂) = 2e, i.e. the row sum 1 of U₁), so H₀ is symmetric
and H₀² = I₄ ⊗ (W₁² + W₂² + 2W₃²). The border entries of W₁² + W₂² + 2W₃² are 4(2n + 1) = 4q and
−e(1 + U₁ ⊗ J₂) + e(1 + U₂ ⊗ J₂) + 2e = 0, and its core block is 4J + 2(U₁² + U₂²) ⊗ J₂ = 2q I ⊗ J₂. So

```
H₀² = 4q I − 2q Λ(I),
```

where for a matrix Y of order 4n, indexed by (k, i), Λ(Y) is the matrix of order 4q with Y[(k, i), (l, j)] D₂[a, b] at
row (k, i, a) and column (l, j, b), and 0 on the borders. (I ⊗ J₂ = 2I − I ⊗ D₂ on the cores.)

**2. Decoupling.** Every row of H₀ is constant on each pair {(l, j, 0), (l, j, 1)}: the borders of the W_k are ±e and
their cores are U ⊗ J₂. Every column of Λ(V) sums to zero over each pair, since the columns of D₂ do. Hence
H₀Λ(V) = 0, and Λ(V)H₀ = (H₀Λ(Vᵀ))ᵀ = 0. Since D₂² = 2D₂, Λ(Y)Λ(Z) = 2Λ(YZ).

**3. The V-part.** Let V be a symmetric (0, ±1) matrix of order 4n with V² = qI whose support is the complement of the
support of the core of H₀: in blocks of order n,

```
P = [[I, J, J, 0],
     [J, 0, I, J],
     [J, I, 0, J],
     [0, J, J, I]]
```

(the core of block (k, l) of H₀ is ±U₁ ⊗ J₂ at P-pattern I, ±U₂ ⊗ J₂ at pattern 0, and 0 at pattern J). Then
H = H₀ + Λ(V) is symmetric, its entries are ±1, and

```
H² = H₀² + H₀Λ(V) + Λ(V)H₀ + Λ(V)² = 4qI − 2qΛ(I) + 0 + 0 + 2Λ(V²) = 4qI.
```

Seberry–Yamada's V-part is the case V = [[I, M, M, 0], [M, 0, −I, −M], [M, −I, 0, M], [0, −M, M, −I]] with M symmetric.

**4. A V from any Hadamard matrix.** Let K be a Hadamard matrix of order n, h = n/2, K = [A | B] with A the first h
columns, E = diag(1, …, 1, −1, …, −1) (h of each), X = [Aᵀ; Aᵀ] and Y = [B, B] (both of order n), and

```
V = [[E,     X,    −EX,   0  ],
     [Xᵀ,    0,     I,    Y  ],
     [−XᵀE,  I,     0,   −YE ],
     [0,     Yᵀ,   −EYᵀ,  E  ]].
```

V is symmetric with support P. Since the columns and the rows of K are orthogonal, AᵀA = BᵀB = nI (of order h),
AᵀB = 0 and AAᵀ + BBᵀ = nI. Hence XXᵀ = YᵀY = [[nI, nI], [nI, nI]] (blocks of order h), which E·E turns into
[[nI, −nI], [−nI, nI]]; XᵀX + YYᵀ = 2(AAᵀ + BBᵀ) = 2nI; XY = 0; and XᵀEX = AAᵀ − AAᵀ = 0, YEYᵀ = 0. The blocks of V²
are then

- (0, 0): E² + XXᵀ + EXXᵀE = I + 2nI = qI, and (3, 3): YᵀY + EYᵀYE + E² = qI;
- (1, 1): XᵀX + I + YYᵀ = qI, and (2, 2): XᵀE²X + I + YE²Yᵀ = qI;
- (0, 1): EX − EX = 0; (0, 2): −E²X + X = 0; (1, 3): −YE + YE = 0; (2, 3): Y − YE² = 0;
- (0, 3): XY + EXYE = 0; (1, 2): −XᵀEX − YEYᵀ = 0.

So V² = qI, and H = H₀ + Λ(V) is a symmetric Hadamard matrix of order 4q. ∎

Remarks.

- Only the support, the symmetry and V² = qI matter. Two Hadamard matrices K₁ = [A₁ | B₁], K₂ = [A₂ | B₂] serve as
  well (X = [A₁ᵀ; A₂ᵀ], Y = [B₁, B₂]), and the blocks I at (1, 2) and (2, 1) may be any diagonal ±1 matrix F (with
  −EXF and −FYE at (0, 2) and (2, 3)).
- Conversely, let V be any symmetric matrix with support P and V² = qI, with diagonal blocks E₁, E₄ at (0, 0), (3, 3)
  and F at (1, 2). Blocks (0, 1) and (2, 3) of V² force V[0, 2] = −E₁XF and V[2, 3] = −FYE₄; then blocks (0, 0) and
  (3, 3) make the rows of X in each sign class of E₁, and the columns of Y in each sign class of E₄, orthogonal; block
  (0, 3) makes those rows orthogonal to those columns when the signs agree; and blocks (1, 1) + (1, 2) give
  X₊ᵀX₊ + Y₊Y₊ᵀ = nI for the rows X₊ of X with E₁ = 1 and the columns Y₊ of Y with E₄ = 1. The columns of
  [X₊ᵀ | Y₊] are orthogonal, so there are at most n of them, and its rank is n, so it is a Hadamard matrix of order n:
  with this U-part, nothing less than a Hadamard matrix of order (q − 1)/2 will do.
- The hypothesis q ≡ 1 (mod 8) is that n = (q − 1)/2 is a multiple of 4; the field is used only for U₁, U₂.
- In [PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md) the propus array needed W₃ = W₄ because the V-part was itself of propus
  form. Here the V-part is not of Williamson type (block (0, 2) of V is −EX ≠ X), so the obstruction there (no
  symmetric orthogonal design of order 4 in four variables) does not apply to it.

## Reach

The census leaves undecided 58 prime powers q ≡ 1 (mod 8) below 2500. [PROPUS_MIYAMOTO.md](PROPUS_MIYAMOTO.md) decided
48 with a symmetric input. The other ten need a Hadamard matrix of order 4m, m = 101, 107, 119, 149, 179, 239, 249, 269,
299 and 305; all are known:

- 404 = 4 · 101, 428 = 4 · 107: Cooper–Wallis with the T-matrices of [t-matrices](../t-matrices) (T-sequences of
  length 101 are Turyn's, 1974; those of length 107 come from Kharaghani and Tayfeh-Rezaie's TT(36), 2005);
- 476 = 4 · 17 · 7, 1196 = 4 · 23 · 13, 1220 = 4 · 5 · 61: Cooper–Wallis with T-matrices of order 17, 23 and 5 (the
  last found by an exhaustive search, ++000, 00+−0, 0000+, 00000) and Turyn's Williamson matrices of order (r + 1)/2
  over GF(r), r = 13, 25 and 121 (I ± B, A, A with Yamada's circulants A, B, A² + B² = rI);
- 596 = 4 · 149, 1076 = 4 · 269: Miyamoto's theorem (1991) with q′ = 149, 269 and Hadamard matrices of orders 148 and
  268 (T-matrices of orders 37 and 67);
- 716 = 4 · 179: Alpöge, Voinov and Reynolds-Haertle (2026; the seeds and border decoded in the JD-Jones-ASES
  Hadamard-B data), a bordered Goethals–Seidel array over ℤ₁₇₈ that verify.py's `bordered-goethals-seidel` builds
  unchanged ([search/data/h716_alpoge.json](../search/data/h716_alpoge.json));
- 956 = 4 · 239: Đoković and Kotsireas's skew matrix ([search/data/skew956_dk.json](../search/data/skew956_dk.json));
- 996 = 4 · 249: [skew/996](../skew/996).

Beyond the census (v ≥ 2500) the theorem gives a symmetric Hadamard matrix of order 4q for every prime power
q ≡ 1 (mod 8) with a Hadamard matrix of order (q − 1)/2; in particular for every such q < 1337, since every multiple
of 4 below 668 is the order of a Hadamard matrix.

## Status

- All ten orders are **New** (to our knowledge, pending the audit of the literature recorded in
  [NOTES.md](../NOTES.md) §30). Their census values are undecided, and the construction of the V-part from a
  non-symmetric Hadamard matrix is not in the sources we read: Seberry–Yamada's Corollary 29 (1990) and its 1992
  restatement (Corollary 8.8) use V = (I, 0, M, M) (or the conference-matrix form) and conclude only ordinary Hadamard
  matrices; Seberry–Balonin's propus paper (2015/2017) applies the propus array only with V = (M, N, I, I); Yamada's
  Theorem 3 (1989), with the same hypothesis as the theorem here, gives ordinary Hadamard matrices through a
  quaternion-type array; Miyamoto's claim of a matrix "of symmetric block form" (1991) needs a symmetric Hadamard
  matrix of order q − 1, and the census says his proof does not give it.
- Ordinary Hadamard matrices of all ten orders are classical (Yamada's Theorem 3, or Miyamoto's theorem with
  H₂ ⊗ K). The symmetry is what is new. Doubling gives symmetric matrices of orders 8q, 16q, ….
- 5732 uses a Hadamard matrix of order 716 that was first constructed in 2026; the other nine inputs were published
  before 2026 except 996 (ours, 2026), for which classical Hadamard matrices of order 996 would serve as well.

The sources searched are listed in [NOVELTY.md](../NOVELTY.md).

## Reproduction and verification

```
python3 search/sym_weighing.py targets          # the census values q = 1 (mod 8) left by sym_propus_miyamoto.py
python3 search/sym_weighing.py certify 809      # writes symmetric/3236/certificate.json
python3 verify/verify.py symmetric/3236         # rebuilds H and checks it
```

`verify/verify.py` (construction `symmetric-propus-weighing`) rebuilds K from its nested certificate with K's own
builder and checks it (H Hᵀ = nI and its SHA-256), builds U₁ and U₂ over GF(q), builds V and checks that it is
symmetric with support P and V² = qI, checks that H₀ vanishes on the support of Λ(V), and checks H = Hᵀ, H Hᵀ = 4qI and
the SHA-256 of the canonical CSV.

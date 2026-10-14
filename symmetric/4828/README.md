# Symmetric Hadamard matrix of order 4828

4828 = 4 · 1207 = 68 · 71, with 1207 = 17 · 71.

## Construction

A symmetric Hadamard matrix is a Hadamard matrix H with H = Hᵀ. This one comes from the symmetric form of a
classical plug-in: a symmetric orthogonal design of type (1, m − 1) and a symmetric Hadamard matrix of order m + 4 give
a symmetric Hadamard matrix of order m(m + 3). Here m = 68. The non-symmetric version — a skew-Hadamard matrix of order
m and the symmetric core of a Hadamard matrix of order m + 4, plugged in with J − 2I — is Wallis's (J. Austral. Math.
Soc. 11 (1970), Corollary 10 (iv); Linear and Multilinear Algebra 3 (1975), Corollary 4, which lists 1207 among new
orders of Williamson-type matrices). Our matrix is a signed row permutation of that one:
P ⊗ X + N ⊗ D = (P ⊗ I)(I ⊗ X + C ⊗ D) with C = PN Paley's skew core. What is new is that the symmetric design makes
it symmetric.

**1. A symmetric OD(68; 1, 67).** Let q = 67, a prime ≡ 3 (mod 4), and let S = I + C be Paley's skew-Hadamard matrix
of order 68 (index ∞ and the elements of GF(67); C[∞, x] = 1, C[x, ∞] = −1, C[x, y] = χ(y − x) with χ the quadratic
character). P = (−1) ⊕ (x ↦ −x), a symmetric signed permutation matrix with P² = I, satisfies PCP = −C (the sign
at ∞ is needed). Hence N = CP is symmetric, has zeros exactly where P is non-zero, N² = 67 I and PN + NP = 0:
xP + yN is a symmetric orthogonal design of type (1, 67) (the same object as in Đoković, arXiv:2508.17141, Theorem 1,
where it is built from difference families).

**2. A core of order 71.** M = H₂ ⊗ M₃₆, with M₃₆ = [[C₁₈ + I, C₁₈ − I], [C₁₈ − I, −C₁₈ − I]] built from Paley's
symmetric conference matrix C₁₈ of order 18 (q = 17), is a symmetric Hadamard matrix of order 72. Moving its first
diagonal entry −1 to position 0 and normalising row and column 0 to (−1, 1, …, 1) by a diagonal ±1 similarity leaves
a symmetric D of order 71 with DJ = J and D² = 72 I − J.

**3. The plug-in.** With X = J − 2I (order 71: −1 on the diagonal, +1 elsewhere),

    H = P ⊗ X + N ⊗ D.

X and D are symmetric and commute (D commutes with J), X² = 4I + 67J and 67 D² = 67 · 72 I − 67 J, so

    H Hᵀ = P² ⊗ X² + (PN + NP) ⊗ XD + N² ⊗ D² = I ⊗ (X² + 67 D²) = (4 + 67 · 72) I = 4828 I,

and H = Hᵀ. Every block of H is ±X or ±D, so H has entries ±1.

**The general statement** (`attack/NOTES.md` §31). If xP + yN is a symmetric OD(m; 1, m − 1) and D is the core of a
symmetric Hadamard matrix of order r + 1, then P ⊗ X + N ⊗ D is a symmetric Hadamard matrix whenever X is a symmetric
±1 matrix commuting with D with X² + (m − 1)D² = m r I. X = J gives r = m − 1 (Đoković's construction, the symmetric
form of Williamson's order k(k − 1)); X = J − 2I gives r = m + 3, since (J − 2I)² = 4I + (r − 4)J must cancel the
−(m − 1)J of (m − 1)D² (the symmetric form of Wallis's m(m + 3)). So symmetric Hadamard matrices of order m(m + 3)
exist whenever a skew-Hadamard matrix of order m with an involutive anti-automorphism (equivalently a symmetric
OD(m; 1, m − 1)) and a symmetric Hadamard matrix of order m + 4 exist. With m = 4n and n odd the order is
4 · n(4n + 3): in the census range this is 1207 (n = 17) and 2185 (n = 23, which needs a symmetric OD(92; 1, 91),
unknown).

## Reproduction

Nothing was searched: the matrix is determined by q = 67 and the symmetric Hadamard matrix of order 72.

```
python3 verify/verify.py symmetric/4828
```

rebuilds it from `certificate.json` (construction `sym-od-plug`: Paley design for q = 67, plug "J-2I", the matrix of
order 72 as Paley II for q = 17 doubled once), checks that xP + yN is a symmetric OD(68; 1, 67), that D is symmetric
with DJ = J and D² = 72 I − J, that X and D commute with X² + 67 D² = 4828 I, and then H Hᵀ = 4828 I, H = Hᵀ and the
SHA-256 of the canonical CSV. The same steps are in `attack/sym_mm1.py` (Paley designs, cores, plug-in) and were
checked on the small cases q = 7, 11, 19, 23 (orders 88, 180, 460, 648).

## Status

To our knowledge this is the first symmetric Hadamard matrix of order 4828. An ordinary (non-symmetric) Hadamard
matrix of this order has been known since Wallis (1970, 1975).

**Novelty checklist** (audit of 2026-10-14 by a research agent; downloaded sources were read as text, never
executed):

1. **The Đoković–Seberry census** (PDF of 2025-05-28, SHA-256 `4b61953c…f04a`): "We list below the 595 odd values of
   v < 2500 for which the existence of symmetric Hadamard matrices of order 4v is undecided: … 1199, 1205, 1207, 1211,
   1213, …". 1207 is in none of its symmetric series (Paley's 2(q + 1) needs q = 2413 = 19 · 127; XXSW needs
   q = 4827 = 3 · 1609; Seberry–Balonin needs q = 1205 or (q + 1)/2 = 603, not prime powers; Turyn C, Turyn, Mathon,
   Seberry–Whiteman and Muzychuk–Xiang fail by arithmetic) and not among its 12 sporadic values. No later version was
   found.
2. **Đoković, arXiv:2508.17141** (the only version): symmetric order q(q + 1) only; no follow-up was found (OpenAlex
   lists no citing work).
3. **The classical non-symmetric constructions** (read in the original papers): Wallis, J. Combin. Theory 6 (1969),
   Corollary 7 (orders (p^k + 1)(p^k + 4)); Wallis, J. Austral. Math. Soc. 11 (1970), Corollary 10 (iv) ("a
   skew-Hadamard matrix of order x − 3 … an Hadamard matrix of order x(x − 3)"); Wallis, Linear and Multilinear
   Algebra 3 (1975), Corollary 4 ("Put X = J − 2I, Y = E … new Williamson type matrices of order 1207, 1827, 2185,
   2575", with E the symmetric core). None of them states a symmetric version, and none of the later surveys or
   papers we read (Seberry–Yamada 1992, Seberry 2013 on strong amicable designs) does.
4. **2026 deposits**: none contains a symmetric matrix of order 4828 (Alpöge et al., Schneider, JD-Jones-ASES,
   Schumer, Kline, Bedani, Karoui, Piza).
5. **Classical routes to a symmetric matrix**: a symmetric conference matrix of order 2414 cannot exist (2413 =
   19 · 127 is not a sum of two squares), which rules out every conference-matrix series; Kronecker products are
   impossible (2414 ≡ 2 (mod 4), 1207 odd); q(q + 1) = 4828 and m(m − 1) = 4828 have no integer solution; regular
   symmetric matrices need a square order.

Not checked: Williamson (1947), the books of Geramita–Seberry (1979), Seberry (2017), Seberry–Yamada (2020) and Hall,
Wallis's Lecture Notes in Mathematics 292, MathSciNet and zbMATH. This is a finding of nothing, not a proof that no
earlier construction exists.

## Files

| File | Contents |
|---|---|
| `certificate.json` | q = 67, the plug-in "J-2I" and the source of the symmetric Hadamard matrix of order 72 (`sym-od-plug`) |
| `H4828.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py symmetric/4828
gzip -dc symmetric/4828/H4828.csv.gz | ./check_matrix
```

SHA-256 of `H4828.csv`: `42a2ec96d4d4630736ecc2837a68f842f8c07324753b9a463150a37f7993af3b`

# Symmetric Hadamard matrix of order 4828

4828 = 4 · 1207 = 68 · 71, with 1207 = 17 · 71.

## Construction

A symmetric Hadamard matrix is a Hadamard matrix H with H = Hᵀ. This one comes from the symmetric form of a classical
plug-in: a symmetric orthogonal design of type (1, m − 1) and a symmetric Hadamard matrix of order m + 4 give a
symmetric Hadamard matrix of order m(m + 3). Here m = 68. The non-symmetric version — a skew-Hadamard matrix of order
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

**The general statement** ([search/NOTES.md](../../search/NOTES.md) §31). If xP + yN is a symmetric OD(m; 1, m − 1)
and D is the core of a symmetric Hadamard matrix of order r + 1, then P ⊗ X + N ⊗ D is a symmetric Hadamard matrix
whenever X is a symmetric ±1 matrix commuting with D with X² + (m − 1)D² = m r I. X = J gives r = m − 1 (Đoković's
construction, the symmetric form of Williamson's order k(k − 1)); X = J − 2I gives r = m + 3, since
(J − 2I)² = 4I + (r − 4)J must cancel the −(m − 1)J of (m − 1)D² (the symmetric form of Wallis's m(m + 3)). So
symmetric Hadamard matrices of order m(m + 3) exist whenever a skew-Hadamard matrix of order m with an involutive
anti-automorphism (equivalently a symmetric OD(m; 1, m − 1)) and a symmetric Hadamard matrix of order m + 4 exist.
With m = 4n and n odd the order is 4 · n(4n + 3). Among the undecided values of the census this gives 1207 (n = 17)
and 2185 (n = 23, which needs a symmetric OD(92; 1, 91), unknown).

## Search

Nothing was searched: the matrix is determined by q = 67 and the symmetric Hadamard matrix of order 72.

## Status

**New** as a symmetric Hadamard matrix. An ordinary Hadamard matrix of order 4828 has been known since Wallis
(1970, 1975); the symmetry is what is new.

- The Đoković–Seberry census lists 1207 among its 595 undecided v for symmetric matrices. 1207 is in none of its
  symmetric series (Paley's 2(q + 1) needs q = 2413 = 19 · 127; XXSW needs q = 4827 = 3 · 1609; Seberry–Balonin needs
  q = 1205 or (q + 1)/2 = 603, not prime powers; Turyn C, Turyn, Mathon, Seberry–Whiteman and Muzychuk–Xiang fail by
  arithmetic) and not among its 12 sporadic values.
- Đoković (arXiv:2508.17141) builds symmetric matrices of order q(q + 1) only. No follow-up was found
  (OpenAlex lists no citing work).
- The non-symmetric constructions are classical: Wallis, J. Combin. Theory 6 (1969), Corollary 7
  (orders (pᵏ + 1)(pᵏ + 4)); Wallis, J. Austral. Math. Soc. 11 (1970), Corollary 10 (iv) ("a skew-Hadamard matrix of
  order x − 3 … an Hadamard matrix of order x(x − 3)"); Wallis, Linear and Multilinear Algebra 3 (1975), Corollary 4
  ("Put X = J − 2I, Y = E … new Williamson type matrices of order 1207, 1827, 2185, 2575", with E the symmetric core).
  None of them states a symmetric version, and neither do the later surveys and papers we read (Seberry–Yamada 1992,
  Seberry 2013 on strong amicable designs).
- No classical route gives a symmetric matrix of order 4828:
  - A symmetric conference matrix of order 2414 cannot exist, since 2413 = 19 · 127 is not a sum of two squares. This
    rules out every conference-matrix series.
  - Kronecker products are impossible (2414 ≡ 2 (mod 4), 1207 odd).
  - q(q + 1) = 4828 and m(m − 1) = 4828 have no integer solution.
  - Regular symmetric matrices need a square order.

We did not read Williamson (1947), Hall's book or Wallis's Lecture Notes in Mathematics 292.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

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

`verify/verify.py` rebuilds the matrix from `certificate.json` (construction `sym-od-plug`: Paley design for q = 67,
plug "J-2I", the matrix of order 72 as Paley II for q = 17 doubled once), checks that xP + yN is a
symmetric OD(68; 1, 67), that D is symmetric with DJ = J and D² = 72 I − J, that X and D commute with
X² + 67 D² = 4828 I, and then H Hᵀ = 4828 I, H = Hᵀ and the SHA-256 of the canonical CSV. The same steps are in
`search/sym_mm1.py` (Paley designs, cores, plug-in) and were checked on the small cases q = 7, 11, 19, 23
(orders 88, 180, 460, 648).

SHA-256 of `H4828.csv`: `42a2ec96d4d4630736ecc2837a68f842f8c07324753b9a463150a37f7993af3b`

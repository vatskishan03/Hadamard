# Symmetric Hadamard matrix of order 5700

5700 = 4 · 1425 = 76 · 75, with 1425 = 3 · 5² · 19.

## Construction

A symmetric Hadamard matrix is a Hadamard matrix H with H = Hᵀ. This one is built like Đoković's symmetric Hadamard
matrices of order q(q + 1) (arXiv:2508.17141, Theorem 2), with q = 75, which is not a prime power, so that theorem
does not apply; the two ingredients are replaced as described below.

**1. A skew-propus quadruple of length 19.** Three ±1 sequences on ℤ₁₉ (−1 at the listed positions)

    A3 = −+++−+−++−−++−+−+++   {0, 4, 6, 9, 10, 13, 15}             symmetric, row sum 5
    A1 = −−−−+++−+++−++−++++   {0, 1, 2, 3, 7, 11, 14}              row sum 5
    A0 = +++++−++−+−+−−+−−−−   {5, 8, 10, 12, 13, 15, 16, 17, 18}   skew type: a0[0] = 1, a0[−g] = −a0[g]

satisfy PAF_A3 + 2 PAF_A1 + PAF_A0 = 76 δ (5² + 2 · 5² + 1² = 76; A0 always has row sum 1). All three are invariant
under the multiplier 7 (of order 3).

**2. A symmetric orthogonal design OD(76; 1, 75).** With A = A3, B = A1, D₀ = A0 − I (circulants) and R the
back-diagonal identity of order 19, put

```
    [ −A     BR      BR      D₀R ]          [ 0   0    0   R ]
N = [  BR    RD₀     A      −RB  ]      P = [ 0   R    0   0 ]
    [  BR    A      −RD₀     RB  ]          [ 0   0   −R   0 ]
    [  D₀R  −RB      RB      A   ]          [ R   0    0   0 ]
```

Then Y = xP + yN is Đoković's matrix Y of arXiv:2508.17141, Theorem 1 (the Balonin array with blocks yA3, yA1, yA1 and
(x − y)I + yA0): P and N are symmetric with complementary supports, P² = I, N² = 75 I and PN + NP = 0, so Y is a
symmetric orthogonal design of type (1, 75). Equivalently, S = I + NP is a skew-Hadamard matrix of order 76 and P an
involutive anti-automorphism of it (PSP = Sᵀ).

**3. The plug-in.** M = P + N is a symmetric Hadamard matrix of order 76. Move its first diagonal entry −1 to position
0 and multiply rows and columns by ±1 so that row and column 0 become (−1, 1, …, 1); the remaining block D (order 75)
is symmetric with DJ = J and D² = 76 I − J. Then

    H = P ⊗ J₇₅ + N ⊗ D

is symmetric (P, N, J and D are) and

    H Hᵀ = P² ⊗ J² + (PN + NP) ⊗ J + N² ⊗ D² = I ⊗ 75 J + 75 I ⊗ (76 I − J) = 5700 I,

because J and D commute. Every block of H is ±J or ±D, so H has entries ±1.

The general statement ([search/NOTES.md](../../search/NOTES.md) §31): a skew-Hadamard matrix of order m with an
involutive anti-automorphism (a symmetric OD(m; 1, m − 1)) gives a symmetric Hadamard matrix of order m(m − 1); D may
be taken from any symmetric Hadamard matrix of order m. Đoković's theorem is the case of the
Xia–Xia–Seberry–Wu difference families for prime powers q = m − 1 ≡ 3 (mod 8), where D is built from the Paley core.
Without symmetry this is Williamson's classical theorem (1944) that a skew-Hadamard matrix of order k gives a Hadamard
matrix of order k(k − 1); an ordinary Hadamard matrix of order 5700 is therefore classical, and only the symmetry is
new.

## Search

`search/sym_mm1.py propus 19 --all` searches skew-propus quadruples of length 19 exhaustively. The row sums must
satisfy a3² + 2 a1² = 75, so (a3, a1) = (5, 5) up to sign. All 84 symmetric A3 and all 512 skew-type A0 are paired,
and the needed autocorrelation of A1 is looked up in a table of all A1 with row sum 5 and a −1 at position 0 (every
such A1 has a shift of this form; 1,311 autocorrelation classes). Up to multipliers, shifts of A1 and A0 ↦ A0ᵀ there
are exactly two solutions, with the same A3 and A0 (7 seconds). The certificate uses the first.

**Reproduction.**

| | |
|---|---|
| Code | commit `2770ccf` (`search/sym_mm1.py`) |
| Command | `python3 search/sym_mm1.py propus 19 --all`; then `python3 search/sym_mm1.py certify -- 19 -+++-+-++--++-+-+++ ----+++-+++-++-++++ +++++-++-+-+--+----` |
| Randomness and caps | none: the search is exhaustive and deterministic |
| Result | two solutions; `sym_mm1.py test` also rebuilds the designs for the Paley primes 3–23 and for every skew-propus length up to 15, and the Szekeres-type designs of `search/szek_pair.c` |

## Status

**New** as a symmetric Hadamard matrix. An ordinary Hadamard matrix of order 5700 is classical (Williamson's theorem,
see Construction); the symmetry is what is new.

- The Đoković–Seberry census lists 1425 among its 595 undecided v for symmetric matrices. 1425 is in none of its eight
  symmetric series nor among its 12 sporadic values.
- Đoković (arXiv:2508.17141) builds symmetric matrices of order q(q + 1) for prime powers q ≡ 3 (mod 8); q = 75 is not
  a prime power. No follow-up was found (arXiv, Crossref, the journal Informatsionno-upravliaiushchie sistemy
  2024–2026), but the forward citations of the paper were not checked in full, since the citation services limited the
  number of requests.
- A symmetric OD(n; 1, n − 1) is the symmetric half of a strong amicable orthogonal design in Seberry's sense
  (Australas. J. Combin. 55 (2013) 5–13, Definition 5). The orders known there are 2ᵗ, pʳ + 1, (n − 1)ʳ + 1 and
  their products, which exclude 76. The blow-up to order n(n − 1) does not appear before Đoković's paper: the 2025
  census still lists 189 and 473 as undecided although Seberry (2013) already has the designs of orders 28 and 44.
- No classical route gives a symmetric matrix of order 5700:
  - A symmetric conference matrix of order 2850 cannot exist, since 2849 = 7 · 11 · 37 is not a sum of two squares.
    This rules out Paley II, Mathon, the Turyn C-series and Seberry–Whiteman at once.
  - XXSW needs q = 5699 = 41 · 139.
  - Seberry–Balonin's Corollary 1 needs q = 1423 ≡ 3 (mod 4).
  - Kronecker products and doubling are impossible for an order ≡ 4 (mod 8) with 2850 ≡ 2 (mod 4).
  - Constant-diagonal and Bush-type matrices need a square order, and two-circulant forms need 5700 to be a sum of two
    squares.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md). The same design also gives amicable Hadamard
matrices of order 76 (S = I + PN and P + N), and through Seberry's 1978 theorem a skew-Hadamard matrix of order 5700:
[skew/5700](../../skew/5700).

## Files

| File | Contents |
|---|---|
| `certificate.json` | the three sequences (positions of −1) and the row sums (`skew-propus-plug`) |
| `sequences.txt` | A3, A1 and A0, one per line |
| `H5700.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py symmetric/5700
gzip -dc symmetric/5700/H5700.csv.gz | ./check_matrix
```

`verify/verify.py` (construction `skew-propus-plug`) rebuilds the sequences, checks the symmetry of A3, the skew type
of A0 and the autocorrelation identity, builds P and N and checks that they form a symmetric OD(76; 1, 75), builds D
and checks DJ = J and D² = 76 I − J, assembles H and checks H Hᵀ = 5700 I, H = Hᵀ and the SHA-256 of the canonical
CSV.

SHA-256 of `H5700.csv`: `feb4a76215d22c517781071defac29c1b01d1f18dc31e6544ce5e64257051ee5`

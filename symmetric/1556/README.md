# Symmetric Hadamard matrix of order 1556

1556 = 4 · 389, with 389 a prime ≡ 5 (mod 8).

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is the case q = 389 of [PERIODIC_GOLAY.md](../PERIODIC_GOLAY.md).
The construction takes Paley's conference matrix of order 390 and a Hadamard matrix K = [[P, Q], [Q, P]] of order 388
from a periodic Golay pair of length 194. H has an involutory automorphism with 776 fixed points.

- **The input.** A periodic Golay pair (a, b) of length 194 over Γ = ℤ₁₉₄ ≅ ℤ₂ × ℤ₉₇ (PAF_a + PAF_b = 388 δ), given in
  `certificate.json` with t = 97, the element of order 2 (the index y of ℤ₁₉₄ is the pair (y mod 2, y mod 97)). It comes
  from our exhaustive search over the pairs invariant under the multiplier group {1, 35, 61} of ℤ₁₉₄* (the subgroup of
  order 3; it acts as the subgroup of order 3 of ℤ₉₇* and trivially on ℤ₂), and is the only pair of this kind up to the
  symmetries listed in [search/pgp/README.md](../../search/pgp/README.md). The row sums are 8 and 18 (8² + 18² = 388),
  with 101 entries +1 in a and 106 in b. K is [[A, B], [−Bᵀ, Aᵀ]] with A, B circulant, reordered by the cosets of ⟨97⟩.
- **The field.** GF(389) = ℤ₃₈₉, χ the quadratic character, ν = 2 (the least non-square).
- **The Q-part.** On the 388 points x ≠ 0: A[x, y] = χ(x + y), Z[x, y] = [x + y = 0],
  F = [[A + Z, A − Z], [A − Z, A + Z]], B = [(𝟙; 𝟙) | (χ; χ) | (k_j; −k_j)] with k_j the columns of K (rows
  relabelled: row r is the point r + 1, row 194 + r the point −(r + 1)), and P = −BᵀFB/1552.
- **The V-part.** V₀[x, y] = χ(xy − 2), V₀[x, ∞] = χ(x), V₀[∞, ∞] = 1, relabelled so that its zeros are the support
  of P.
- **The matrix.** H = [[F, B, B], [Bᵀ, P + V, P − V], [Bᵀ, P − V, P + V]].

## Search

The periodic Golay pair of length 194 comes from an exhaustive enumeration of the pairs (a, b) of ±1 sequences on
ℤ₂ × ℤ₉₇ that are constant on the orbits of H = {1, 35, 61} (33 free signs per layer and sequence, 132 in all). With
A = (a₀ + a₁)/2 and C = (a₀ − a₁)/2 for the two layers of a, and B, D likewise from b, the pair is a periodic Golay pair
exactly when AA* + BB* = 97 δ and CC* + DD* = 97 δ in ℤ[ℤ₉₇], where A, C (and B, D) are ternary with complementary
supports. Stage 1 lists the rows (X, Y), pairs of H-invariant ternary sequences with XX* + YY* = 97 δ, up to the
symmetries; stage 2 joins the rows whose supports are complementary. The sums of X and Y over ℤ₉₇ are forced to 4 and 9
(97 = 4² + 9², and a sum is congruent to the value at 0 modulo 3). Stage 1 enumerated 8,894,476,872,902 sign patterns in
24,247 buckets and found 897,576 rows; stage 2 found 2 matches, which are one equivalence class of pairs, the pair of
`certificate.json`.

The run used 8 threads of a 32-core AMD EPYC 9R14 (246 GB) in four sessions, resumed from a checkpoint after three
interruptions: 46,419 s (12.9 h) of wall time in all by the scheduler's count (the four sessions add up to 46,430 s), 99.4 thread-CPU hours, about 1.1 GB of memory, from 1 to 2
November 2026. The search is exhaustive for this family only: it covers no pair with another multiplier group or none.
The programs, the exact command and the validation are in [search/pgp](../../search/pgp/README.md) and
[search/NOTES.md](../../search/NOTES.md) §113.

The pair is stored in `certificate.json` and hard-coded in `search/sym_involution.py` (`PAIRS[389]`), and
`python3 search/sym_involution.py certify 389` builds the certificate from it. The certificate alone proves the result.

## Status

**New** as a symmetric Hadamard matrix.

- The Đoković–Seberry census (2025-05-28) lists 389 among its 595 undecided v for symmetric matrices, in none of its
  series (Paley needs q = 777 = 3 · 7 · 37, Xia–Xia–Seberry–Wu q = 1555 = 5 · 311, Seberry–Balonin q = 387 = 3² · 43) and
  not among its sporadic values. [PROPUS_WEIGHING.md](../PROPUS_WEIGHING.md) needs q ≡ 1 (mod 8) and does not apply to
  q = 389. Đoković's orders q(q + 1) (arXiv:2508.17141) do not reach it: q(q + 1) = 1556 has no integer solution.
  Handbook Table V.1.52 stops at v < 200, so the census is the reference list here.
- The order is within the scope of Miyamoto's claim, as 404 and 596 are. Miyamoto (J. Combin. Theory Ser. A 57 (1991),
  Theorem 5, last sentence, as reported by the census) claimed a Hadamard matrix of order 4q "of symmetric block form"
  whenever a symmetric conference matrix of order q + 1 and a symmetric Hadamard matrix of order q − 1 exist, and Seberry
  and Yamada restate it (JCMCC 7 (1990), Corollary 35; the 1992 survey, Corollary 8.11). For q = 389 both inputs exist:
  Paley's symmetric conference matrix of order 390 and a symmetric Hadamard matrix of order 388 (v = 97 is in the
  census's Paley series: Paley II with q = 193). The census says the claim is ignored in Miyamoto's proof, that the
  diagonal blocks he displays are not symmetric, and it lists 389 as undecided
  ([PROPUS_MIYAMOTO.md](../PROPUS_MIYAMOTO.md#status)). We could not read his paper and do not use the claim; the
  construction here has its own proof ([PERIODIC_GOLAY.md](../PERIODIC_GOLAY.md)) and a different input.
- Ordinary Hadamard matrices of order 1556 are classical: Miyamoto's theorem with Paley II of order 388 (the Hadamard
  Matrix Explorer builds 1556 as `M(P2(388))`), and, as far as a damaged text extraction of Seberry–Yamada's Table A.1
  (1992) can be read, Williamson-type matrices of order 389 (entry 389, key wk). The symmetry is what is new.
- The input is new too. Periodic Golay pairs of length 194 were undecided in Đoković and Kotsireas's list of lengths
  below 300 (Numer. Algorithms 69 (2015) 523–530, arXiv:1310.5773, Section 5), and the smallest undecided length was
  still 106 in Lumsden, Kotsireas and Bright (arXiv:2408.15611, Math. Comp. 95 (2026)). We found no source that gives a
  pair of length 194; the pair here is, to our knowledge, the first. It comes from one multiplier group, so the search
  does not decide the existence of other pairs of length 194 ([search/NOTES.md](../../search/NOTES.md) §113).
- Not excluded: that a known Hadamard matrix of order 388 (Paley II, Đoković's skew matrices of order 388) already has
  the form K = [[P, Q], [Q, P]], which would give this matrix without the pair. We did not test it.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, the group, t, the periodic Golay pair (a, b), and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/1556
```

`verify/verify.py` (construction `symmetric-involution`) checks the pair and K, builds the Q-part and V and checks
their identities, assembles H and checks H = Hᵀ, H Hᵀ = 1556 I and the SHA-256 of the canonical CSV.

SHA-256 of `H1556.csv`: `3eea11362d17330271f31be78129dc596364c4ba99a6540f01abe7584e3c85ec`

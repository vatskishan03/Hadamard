# Hadamard matrix of order 4204

4204 = 4 · 1051, with 1051 prime.

## Construction

The same three steps as [6796](../6796), with the notation defined there and the bordered form of the core. Each step
is checked exactly by `verify/verify.py`.

**1. A quaternary conference matrix of order 526.** The bordered zero-diagonal Szöllősi array
([search/NOTES.md](../../search/NOTES.md) §75) over ℤ₂₆₂ = ℤ₂ × ℤ₁₃₁:
- A and B map ℤ₂₆₂ to {0, ±1}, are symmetric and vanish exactly at 0; C and D map it to {±1};
- PAF_A(g) + PAF_B(g) + PAF_C(g) + PAF_D(g) = −4 for g ≠ 0;
- the row sums are 1, 1, 0, 0.

With ω = (1 + i)/2, the interior blocks are X′ = ω X_A + ω̄ X_B and Y′ = ω X_C + ω̄ X_D R, with R[g, −g] = 1. Here
ρ = ω r_A + ω̄ r_B = 1, and the border units are t = v = 1. The blocks are

```
X = [ 0  eᵀ ]      Y = [ −1  eᵀ ]
    [ e  X′ ]          [  e  Y′ ]
```

and C = [[X, Y], [Y*, −X*]] has order 526, zero diagonal, entries ±1, ±i elsewhere, and C C* = 525 I.

The sequences are constant on the orbits of H = ⟨53⟩, of order 5: 53 is odd and has order 5 modulo 131. C and D take
one value on each of the 54 orbits of H, and A and B on each of the 28 orbits of ⟨H, −1⟩.

**2. A quaternary Hadamard matrix of order 2102.** q = 2 · 526 − 3 = 1049 is a prime ≡ 1 (mod 8), and the complex
Yamada–Kiyasu template ([NOTES.md](../../NOTES.md) §31, §34) gives K of order 2102 with K K* = 2102 I
([complex-hadamard/C2102.json](../../complex-hadamard/C2102.json)).

**3. Turyn's doubling** gives order 4204.

## Search

This structure had the highest estimate among those examined for the open orders of this family (10^8.9 solutions,
calibrated on analogues of order 2p), but it could not be searched with lists of full sequences: the pools for C and D
hold about 10¹⁴ sequences.

The matrix was found by a two-level sampler (compress, sample the compressed instances, complete each one exactly
through its fibres), as in the FLC classification method ([17] in the [README](../../README.md#references)). The sums
c = x(0, ·) + x(1, ·) and differences d = x(0, ·) − x(1, ·) over the ℤ₂ factor split the condition exactly:

- **Level 1:** Σ PAF_c = 1050 δ − 8 on ℤ₁₃₁ (1042 at 0, −8 elsewhere). Away from 0, c takes values in {−2, 0, 2}.
- **Level 2:** Σ PAF_d = 1050 δ. Away from 0, d is 0 exactly where c is ±2, so d only fixes signs on the zeros of c.

The level-1 instances were sampled in the stratum z_C = z_D = 18, the number of zero orbits of C and D. The binomial
law of the full solutions favours that stratum over the uniform population. Each instance was completed exactly. The
first solution came 103 seconds after launch; 19 were found, 17 of them distinct
([search/NOTES.md](../../search/NOTES.md) §78).

The certificate alone proves the result; the search programs for this order are not included in this repository.

## Status

**New.** The entry 1051(3) is unresolved in Cati–Pasechnik Table 4 (v1 and v2; Handbook Table 1.53 also gives t = 3),
and Table 5 of Egan, Kotsireas, London and Ó Catháin (arXiv:2610.27175, 19 October 2026) lists 1051 as open. No later
construction was found. Table 4 already knows 8 · 1051, so doubling settles nothing further. Đoković's set Δ and his
conversions of 2010 do not contain 1051, and SageMath's development version (20 October 2026) has no branch that
reaches 4204. None of the standard recursive constructions applies:

- Paley: 4203 = 3² · 467 and 2101 = 11 · 191 are not prime powers. No symmetric conference matrix of order 2102
  exists, since 2101 is not a sum of two squares.
- Miyamoto: 1051 ≡ 3 (mod 4).
- Yamada: it needs a symmetric conference matrix of order 526, which cannot exist, because 525 = 3 · 5² · 7 is not a
  sum of two squares.
- Kiyasu: q = 1049 ≡ 1 (mod 8), not 5.
- Product theorems: 4204 ≡ 4 (mod 8).
- Cooper–Wallis: 1051 is prime, so it needs T-matrices or Williamson-type matrices of order 1051. Neither is known.
- The block routes of [PALEY_BLOCKS.md](../PALEY_BLOCKS.md) and [MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md): they need
  4203 a prime power or 1051 ≡ 1 (mod 4). Paley's matrix of order 1052 serves n = 263 · t there, not 1051.

**Credits.** As for [3436](../3436/README.md#status). The zero-diagonal arrays, plain and bordered, are ours
([search/NOTES.md](../../search/NOTES.md) §74–75). Miyamoto's Corollary 4 at m = 525 would need a quaternary
conference matrix of order 526 with Cᵀ = C; the core here has Cᵀ ≠ C as built.

**Ingredients.**
- The quaternary conference matrix of order 526 is new to our knowledge. No real one exists: 525 is not a sum of two
  squares, and 526 ≡ 2 (mod 4). The families of complex ones need 525 to be a prime power. Fallon and Iverson's
  ancillary library (arXiv:2610.09169) has a certified conference signal over ℤ₅₂₅, which gives a complex conference
  matrix of order 526. Its invariants C_ij C_kl C̄_il C̄_kj reach 0.77 from {±1, ±i}, so it is not equivalent to a
  quaternary one.
- The quaternary Hadamard matrix of order 2102 lies above Seberry–Yamada's table (N < 1000) and Zaitsev et al.'s page
  (orders below 2000). Any such matrix gives H(4204) by Turyn's doubling, so none was known while 4204 was open. Its
  other routes fail: 4203 and 2101 are not prime powers, no Williamson or Williamson-type matrices of order 1051 are
  known, and ℓ = 1050 lies in neither of Jedwab and Pender's families of quaternary Legendre pairs.

**Risk.** Zaitsev et al.'s page stops below 2000. Middletownbooks/Hadamard does not reach 1051; its open list stops at
971, and its last commit (20 October) changed only the README. The Seberry–Yamada book of 2020 was not read; by Cati
and Pasechnik's account its Table A.17 also has m > 2 at 1051. Miyamoto's lists of 1991 were not read. The sources
searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | the group, the border, the symmetry group, the four sequences, and q, r, ξ for Yamada's circulants |

The matrix is not stored. `python3 verify/verify.py --write 4204` writes `H4204.csv`.

## Verification

```
python3 verify/verify.py 4204
python3 verify/verify_complex.py 2102
```

SHA-256 of `H4204.csv`: `35b69b50dd0999be8287dbde1978f602aa976a25c7d9ba4b755fa6f2a40c4c1d`

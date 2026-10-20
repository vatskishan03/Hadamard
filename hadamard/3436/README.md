# Hadamard matrix of order 3436

3436 = 4 · 859, with 859 prime.

## Construction

The same three steps as [6796](../6796), each checked exactly by `verify/verify.py`.

**1. A quaternary conference matrix of order 430.** The zero-diagonal Szöllősi array
([search/NOTES.md](../../search/NOTES.md) §74–76) over the cyclic group ℤ₂₁₅ = ℤ₅ × ℤ₄₃. A and B map ℤ₂₁₅ to
{0, ±1}, are symmetric and vanish exactly at 0; C and D map it to {±1}; and

```
PAF_A(g) + PAF_B(g) + PAF_C(g) + PAF_D(g) = 0   for g ≠ 0.
```

With ω = (1 + i)/2, X_f[g, h] = f(h − g) and the reversal R[g, 214 − g] = 1, the matrices X = ω X_A + ω̄ X_B and
Y = ω X_C + ω̄ X_D R give C = [[X, Y], [Y*, −X*]] of order 430 with zero diagonal, entries ±1, ±i elsewhere, and
C C* = 429 I.

The sequences are constant on the orbits of H = ⟨49⟩, of order 6: 49 ≡ −1 (mod 5) and 49 ≡ 6 (mod 43), where 6 has
order 3. Every block is therefore symmetric in the ℤ₅ coordinate. The row sums are 14, 18, 13 and 13:
14² + 18² + 13² + 13² = 858 = 4 · 215 − 2.

**2. A quaternary Hadamard matrix of order 1718.** q = 2 · 430 − 3 = 857 is a prime ≡ 1 (mod 8), and the complex
Yamada–Kiyasu template ([NOTES.md](../../NOTES.md) §31, §34) gives K of order 2(q + 2) = 1718 with K K* = 1718 I
([complex-hadamard/C1718.json](../../complex-hadamard/C1718.json)).

**3. Turyn's doubling** H = [[K_r + K_i, K_r − K_i], [−(K_r − K_i), K_r + K_i]] gives order 3436.

## Search

The structure (ℤ₂₁₅, ⟨49⟩) was the best-ranked one for this order in rounds 20 and 21 (search/NOTES.md §75). A
meet-in-the-middle join of full sequences over its 200,000 most likely split values found nothing.

This matrix was found by a two-level sampler in the manner of FLC ([github.com/vatskishan03/FLC](https://github.com/vatskishan03/FLC)):
- **Level 1** is the compression to ℤ₄₃: per column, g₀ = f₀ + 2f₁ + 2f₂. The level-1 instances are compressed
  quadruples with Σ PAF = 858 δ on ℤ₄₃, conditioned on the ℤ₅ row tuple.
- **Level 2** chooses the signs of v = (f₁ − f₂)/2 on the columns with |g₀| = 1. Each key group of level-1 pairs is
  completed at once by a hashed join.

The first solution came after about two CPU-minutes, from 1/64 of each of the large level-1 lists. Four distinct
cores were found, in two row patterns (search/NOTES.md §78).

## Status

**New.** The entry 859(3) is unresolved in Cati–Pasechnik Table 4 (v1 and v2; Handbook Table 1.53 also gives t = 3),
and no later construction was found. Table 4 already knows 8 · 859, so doubling settles nothing further. Đoković's set
Δ and his conversions of 2010 do not contain 859. SageMath's development version (20 October 2026) has no branch that
reaches 3436: the only branch not closed by arithmetic needs an srg(1717, 858, 428, 429), which no family in Sage
builds and which lies above Brouwer's table. None of the standard recursive constructions applies:

- Paley: 3435 = 3 · 5 · 229 and 1717 = 17 · 101 are not prime powers, and no symmetric conference matrix of order 1718
  is known.
- Miyamoto: 859 ≡ 3 (mod 4).
- Yamada: it needs a symmetric conference matrix of order 430, which cannot exist, because 429 = 3 · 11 · 13 is not a
  sum of two squares.
- Kiyasu: q = 857 ≡ 1 (mod 8), not 5.
- Product theorems: 3436 ≡ 4 (mod 8).
- Cooper–Wallis: 859 is prime, so it needs T-matrices or Williamson-type matrices of order 859. Neither is known.
- The block routes of [PALEY_BLOCKS.md](../PALEY_BLOCKS.md) and [MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md): they need
  3435 a prime or 859 ≡ 1 (mod 4).

**Credits.** As for [6796](../6796/README.md#status). Miyamoto's Corollary 4 (1991, as transcribed in the
JD-Jones-ASES Hadamard-M note) reaches the same orders 4(q + 2) from a C₂-matrix of order q + 3. In his displayed
form, such a matrix is the same thing as a quaternary conference matrix C of order (q + 3)/2 with Cᵀ = C. The template
needs no such symmetry, and the core here has Cᵀ ≠ C as built. The sampler adapts FLC (fibre-lattice completion,
[github.com/vatskishan03/FLC](https://github.com/vatskishan03/FLC)), an exhaustive classifier, to sampling.

**Ingredients.**
- The quaternary conference matrix of order 430 is new to our knowledge. No real one exists: a symmetric one would
  need 429 to be a sum of two squares, and a skew one needs order ≡ 0 (mod 4). The families of complex ones
  (Paley-type, Berman, BGW and GC, Seberry–Whiteman, Et-Taoui) need 429 to be a prime power. Fallon and Iverson's
  ancillary library (arXiv:2610.09169) has a certified conference signal over ℤ₄₂₉, which gives a complex conference
  matrix of order 430. That matrix is not equivalent to a quaternary one: some of its invariants C_ij C_kl C̄_il C̄_kj
  lie 0.77 from {±1, ±i}, the largest distance possible.
- The quaternary Hadamard matrix of order 1718 settles an order that is blank in Seberry–Yamada's Table 11.2 (1992)
  and on the unknown list of Zaitsev et al.'s page (last updated 14 October 2026, unchanged on 20 October). Its other
  routes fail: 3435 and 1717 are not prime powers (Ito–Egan, Turyn), no Williamson or Williamson-type matrices of
  order 859 are known, and ℓ = 858 lies in neither of Jedwab and Pender's families of quaternary Legendre pairs.

**Risk.** Zaitsev et al. use the same template and do not date their matrices; their page stops below 2000 and lists
1718 as unknown. Middletownbooks/Hadamard lists 859 as work in progress: a Goethals–Seidel structure over ℤ₈₅₉ with a
multiplier group of order 13, judged out of reach. Its last commit is of 18 October. The Seberry–Yamada book of 2020
was not read; by Cati and Pasechnik's account its Table A.17 also has m > 2 at 859. Miyamoto's lists of 1991 were not
read; his Corollary 4 at m = 429 would need a symmetric quaternary conference matrix of order 430, which is not known.
The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | the group, the reversal, the symmetry group, the four sequences, and q, r, ξ for Yamada's circulants |

The matrix is not stored. `python3 verify/verify.py --write 3436` writes `H3436.csv`.

## Verification

```
python3 verify/verify.py 3436
python3 verify/verify_complex.py 1718
```

SHA-256 of `H3436.csv`: `752c0014cffb470a2bb1a83e0ab4f3842b9d96785b698dd8f3a1778331b43f70`

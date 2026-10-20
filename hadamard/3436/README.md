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

**New.** The entry 859(3) is unresolved in Cati–Pasechnik Table 4, and no later construction was found. Table 4
already knows 8 · 859, so doubling settles nothing further. None of the standard recursive constructions applies:

- Paley: 3435 = 3 · 5 · 229 and 1717 = 17 · 101 are not prime powers.
- Miyamoto: 859 ≡ 3 (mod 4).
- Yamada: it needs a symmetric conference matrix of order 430, which cannot exist, because 429 = 3 · 11 · 13 is not a
  sum of two squares.
- Kiyasu: q = 857 ≡ 1 (mod 8), not 5.
- Product theorems: 3436 ≡ 4 (mod 8).

The quaternary Hadamard matrix of order 1718 settles an order that Zaitsev et al.'s page (last updated 14 October
2026) lists as unknown. Credits for the array and the template are as for [6796](../6796/README.md#status). The
sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

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

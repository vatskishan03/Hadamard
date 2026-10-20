# Hadamard matrix of order 4204

4204 = 4 · 1051, with 1051 prime.

## Construction

The same three steps as [6796](../6796), with the bordered form of the core. Each step is checked exactly by
`verify/verify.py`.

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

Round 21 rated this structure the strongest live one for any of the targets (estimate 10^8.9 solutions, calibrated on
analogues of order 2p), but it could not be searched with lists of full sequences: the pools for C and D hold
about 10¹⁴ sequences.

The matrix was found by a two-level sampler in the manner of FLC
([github.com/vatskishan03/FLC](https://github.com/vatskishan03/FLC)). The sums c = x(0, ·) + x(1, ·) and differences
d = x(0, ·) − x(1, ·) over the ℤ₂ factor split the condition exactly:

- **Level 1:** Σ PAF_c = 1042 δ − 8 on ℤ₁₃₁. The values of c lie in {−2, 0, 2}.
- **Level 2:** Σ PAF_d = 1050 δ. Here d is 0 exactly where c is ±2, so d only fixes signs on the zeros of c.

The level-1 instances were sampled in the stratum z_C = z_D = 18, the number of zero orbits of C and D. The
binomial law of the full solutions favours that stratum over the uniform population. Each instance was completed
exactly. The first solution came 103 seconds after launch; 19 were found, 17 of them distinct
(search/NOTES.md §78).

## Status

**New.** The entry 1051(3) is unresolved in Cati–Pasechnik Table 4, and no later construction was found. Table 4
already knows 8 · 1051, so doubling settles nothing further. None of the standard recursive constructions applies:

- Paley: 4203 = 3² · 467 and 2101 = 11 · 191 are not prime powers.
- Miyamoto: 1051 ≡ 3 (mod 4).
- Yamada: it needs a symmetric conference matrix of order 526, which cannot exist, because 525 = 3 · 5² · 7 is not a
  sum of two squares.
- Kiyasu: q = 1049 ≡ 1 (mod 8), not 5.
- Product theorems: 4204 ≡ 4 (mod 8).

The quaternary Hadamard matrix of order 2102 lies above the range of the 1992 table and of Zaitsev et al.'s page. The
bordered zero-diagonal array is ours (search/NOTES.md §75); credits for the template are as for
[6796](../6796/README.md#status). The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

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

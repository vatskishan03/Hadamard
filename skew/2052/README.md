# Skew-Hadamard matrix of order 2052

2052 = 4 · 513, with 513 = 3³ · 19.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Group.** The difference family lives on G = ℤ₁₉ × GF(27)⁺ = ℤ₁₉ × ℤ₃ × ℤ₃ × ℤ₃, a non-cyclic group of order 513.
Elements are quadruples (a, b₀, b₁, b₂) in mixed-radix order (last coordinate fastest); (b₀, b₁, b₂) are the
coordinates of b = b₀ + b₁x + b₂x² in GF(27) = F₃[x]/(f), with f the first monic irreducible cubic in the search order
of `gf_model` in `search/paley_config.py`.

**Symmetry: a Paley product.** All four sequences are invariant under H = the 13 squares of GF(27)*, acting by
multiplication on the GF(27) factor and trivially on ℤ₁₉. On coordinates its generator
(multiplication by g², g primitive) is the matrix stored in `certificate.json` under `automorphisms`. Since
27 ≡ 3 (mod 4), −1 is not a square, so x ↦ −x pairs the orbits and a skew first block is possible.

Every block is u ⊗ 1 + w ⊗ J + z ⊗ χ, where χ is the quadratic character of GF(27), J the all-ones function on GF(27),
and u, w, z are sequences on ℤ₁₉ (see [search/NOTES.md](../../search/NOTES.md) §15).

**Difference family.** x₁, …, x₄ are ±1 functions on G whose periodic autocorrelations sum to zero at every nonzero
shift. The first is of skew type: x₁(0) = 1 and x₁(−g) = −x₁(g). The row sums are 1, 27, 19, 31, and
1² + 27² + 19² + 31² = 2052.

**The matrix.** With group-developed blocks X[g, h] = x(h − g) and R[g, h] = 1 when g + h = 0, the Goethals–Seidel
array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is a Hadamard matrix of order 2052 with H + Hᵀ = 2I.

## Search

1. **Reduction.** The compressions y_k = u_k + 27 w_k of the blocks to ℤ₁₉ satisfy Σ_k Σ_a y_k(a)² = 4 · 19 · 27,
   which reads 28 N_w − 4 N₋ = 76: three positions carry w ≠ 0 (N_w = 3, N₋ = 2). Two of them in one block would leave
   a term of size 27² in the autocorrelations of the compressions, so they sit one in each of the blocks B, C, D.
2. **Exact two-stage search** (`search/paley_config.py --field 3,3`). Stage 1 lists all 12,124,320 quadruples of
   compressions to ℤ₁₉ whose autocorrelations sum to zero. Stage 2 takes the first and joins all z-quadruples with the
   prescribed autocorrelation sums and zero cross terms (`search/join4x.c`): 6 solutions. The first two were lifted to
   ℤ₁₉ × GF(27) and verified.

**Reproduction.**

| | |
|---|---|
| Code | commit `520f0b6` |
| Command | `python3 search/paley_config.py 513 --field 3,3 --maxu 6` |
| Randomness and caps | none in the search (exact enumeration); which solutions join4x prints first depends on the OpenMP schedule |
| Result | 12,124,320 y-quadruples (stage 1: 36 minutes); the first extends to 6 z-solutions (stage 2: about 30 minutes) |

## Status

**New** as a skew-Hadamard matrix, and so is its doubling, of order 4104 = 8 · 513. An ordinary Hadamard matrix of
order 2052 is known (Cooper–Wallis, listed by Piza on Zenodo as "2052 = CW (19, 27)"). It is not skew.

- Cati–Pasechnik Table 5 lists 513(4) in plain type, an entry inherited from Seberry–Yamada (2020, Table 9.2) and not
  updated. Seberry–Yamada (1992, Table 7.1) has "513 4". The smallest known skew order of the form 2ᵐ · 513 was 8208,
  so 2052 and 4104 were unknown.
- The Đoković–Seberry census lists 513 among its 737 undecided v. 513 is in none of its Paley, Szekeres,
  Szekeres–Whiteman or Spence series nor among its 41 sporadic values.
- SageMath returns Unknown for 2052. It stores no difference family on any group of order 513.
- No infinite family reaches it:
  - Paley I (2051 = 7 · 293 is not a prime power).
  - Szekeres, Szekeres–Whiteman, Spence and Momihara–Xiang 2(qᵉ + 1) (qᵉ = 1025 = 5² · 41 is not a prime power).
  - Wallis–Whiteman and Momihara–Xiang 4(qᵉ + 1) (qᵉ = 512 is even).
  - Spence 4(q² + q + 1) (2049 is not a square).
  - (p − 1)ᵘ + 1.
  - Products 2ᵗ ∏(qᵢ + 1), skew × amicable and other product theorems (2052 has only two factors of 2).
- One classical theorem has admissible parameters: Seberry (Ars Combin. 6 (1978), Thm 7) with n = 4, m = 171 would
  give 2052 from an OD(684; 1, 171, 512). We found no such design recorded.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | the group ℤ₁₉ × ℤ₃³, the automorphisms generating H, the row sums and the four sequences (`goethals-seidel-group`) |
| `sequences.txt` | the four sequences, one per line, elements in mixed-radix order |
| `H2052.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/2052
gzip -dc skew/2052/H2052.csv.gz | ./check_matrix
```

`verify/verify.py` also checks the invariance under the stated automorphism.

SHA-256 of `H2052.csv`: `59e1f17a54aaee2d0ceaf6d94dc46385c2bacb75519145d4548827afdaface40`

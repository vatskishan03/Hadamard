# Skew-Hadamard matrix of order 1372

1372 = 4 · 343, with 343 = 7³.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Group.** The difference family lives on the additive group of the field GF(343) = GF(7³), that is
G = ℤ₇ × ℤ₇ × ℤ₇. Elements are triples in mixed-radix order (last coordinate fastest). The cyclic group
ℤ₃₄₃ has no comparable structure.

**Symmetry: cyclotomy.** Let ω be a primitive element of GF(343), with multiplicative order 342. Multiplication by
ω is a linear map of ℤ₇³, the companion matrix of a primitive cubic. All four sequences are invariant
under H = ⟨ω¹⁸⟩, the multiplicative subgroup of order 19. On coordinates its generator is the matrix
[[2, 1, 6], [0, 6, 4], [3, 4, 2]] (mod 7). So every block is a union of cyclotomic classes of order 18, together
with the point 0.

H has odd order, so x ↦ −x pairs the 18 classes, and a skew first block is possible. There are 19 orbits
and 9 autocorrelation classes.

**Difference family.** x₁, …, x₄ are ±1 functions on GF(343) whose periodic autocorrelations sum to zero
at every nonzero shift. The first is of skew type: x₁(0) = 1 and x₁(−g) = −x₁(g). The row sums are
1, 37, 1, 1, and 1² + 37² + 1² + 1² = 1372.

**The matrix.** With group-developed blocks X[g, h] = x(h − g) and R[g, h] = 1 when g + h = 0, the
Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is a Hadamard matrix of order 1372 with H + Hᵀ = 2I.

## Search

1. **Rank.** `search/skew_field.py rank 7 3` enumerates the odd-order cyclotomic subgroups of GF(343)*.
   For order 19 the Gaussian model predicts about 3·10⁵ solutions:
   - 512 skew first blocks;
   - pools of about 5·10⁴ sequences for the other three blocks.
2. **Match.** `search/skew_field.py run 7 3 19 1,1,37` matches the four lists exactly (`find4b`, pair
   work 2·10⁹).
3. **Result.** It found 4,211,424 solutions in a few minutes.

**Reproduction.**

| | |
|---|---|
| Code | commit `f904af6` (the options used behave the same at later commits) |
| Command | `python3 search/skew_field.py run 7 3 19 1,1,37` |
| Randomness and caps | exact (no sampling, no seed); pair-work cap 6·10¹⁰ |
| Result | 4,211,424 solutions |

## Status

**New** as a skew-Hadamard matrix. Cati–Pasechnik Table 5 lists 343(4): the smallest known skew order of the form
2ᵐ · 343 was 2⁴ · 343 = 5488, so 1372 and 2744 were unknown. The order is above 1208, the range Cati and Pasechnik
checked against implemented constructions.

- Đoković's papers of 1992–2024 (alone and with Kotsireas, Golubitsky or Balonin), Koukouvinos–Stylianou (Discrete
  Math. 2008, Tables 1 and 3) and the Đoković–Seberry census do not construct it. The census lists 343 among its 737
  undecided values.
- SageMath's `unknown_skew_hadamard` list stops at 1196, below this order.
- The known families do not reach 1372:
  - the skew Hadamard difference families of Momihara and Xiang need q ≡ 2ᵘ + 1 (mod 2ᵘ⁺¹) with u ≥ 2, but 343 ≡ 7
    (mod 16);
  - the Szekeres construction gives order 2(q + 1);
  - the product theorem would need 1372 = n(pᵏ + 1) with a skew order n;
  - Paley I, Wallis–Whiteman, Spence and amicable products were also checked.

The abundance of solutions (more than 4·10⁶ in this one structure) suggests that cyclotomic Goethals–Seidel families
over GF(p³) and GF(p²) have been searched less than cyclic ones.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md). Doubling (Seberry) gives a skew-Hadamard matrix of
order 2744.

## Files

| File | Contents |
|---|---|
| `certificate.json` | the four sequences as ± strings over ℤ₇³, the generator of H as a matrix over ℤ₇, the row sums |
| `sequences.txt` | the four sequences, one per line |
| `H1372.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/1372
gzip -dc skew/1372/H1372.csv.gz | ./check_matrix
```

SHA-256 of `H1372.csv`: `32401cf85e1d4b1284c1def607db2fa1e89e90d82babde46075c73478273f5c5`

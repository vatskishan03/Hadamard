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

1. **Rank.** `attack/skew_field.py rank 7 3` enumerates the odd-order cyclotomic subgroups of GF(343)*.
   For order 19 the Gaussian model predicts about 3·10⁵ solutions:
   - 512 skew first blocks;
   - pools of about 5·10⁴ sequences for the other three blocks.
2. **Match.** `attack/skew_field.py run 7 3 19 1,1,37` matches the four lists exactly (`find4b`, pair
   work 2·10⁹).
3. **Result.** It found 4,211,424 solutions in a few minutes.

## Reproduction

| | |
|---|---|
| Code | commit `58b1460` of this branch (the options used behave the same at later commits) |
| Command | `python3 attack/skew_field.py run 7 3 19 1,1,37` |
| Randomness and caps | exact (no sampling, no seed); pair-work cap 6·10¹⁰ |
| Result of that run | 4,211,424 solutions |

The tools write their hit logs, temporary list files and compiled helpers to `$HADAMARD_WORK` (default
`<tmp>/hadamard-work`), never into the repository. `attack/skew_certify.py` (or `skew_group.py certify` for the
non-cyclic orders) turns a logged hit into the files in this directory. The certificate does not depend on
the search: `verify/verify.py` rebuilds and checks the matrix from `certificate.json` alone.

## Status

To our knowledge this is the first skew-Hadamard matrix of order 1372.

- Cati and Pasechnik, *A database of constructions of Hadamard matrices*, arXiv:2411.18897v2 (Aug. 2025), Table 5, lists 343(4). So the smallest known skew order of the form 2ᵐ · 343 was 5488, and 1372 and 2744 were unknown.
- Doubling now gives 2744 as well.
- The known families do not reach 1372:
  - the skew Hadamard difference families of Momihara and Xiang need q ≡ 2ᵘ + 1 (mod 2ᵘ⁺¹) with u ≥ 2, but 343 ≡ 7 (mod 16);
  - the Szekeres construction gives order 2(q + 1);
  - the product theorem would need 1372 = n(pᵏ + 1) with a skew order n.
- A literature search in October 2026 found no construction.

The abundance of solutions (more than 4·10⁶ in this one structure) suggests that cyclotomic
Goethals–Seidel families over GF(p³) and GF(p²) have been searched less than cyclic ones.

## Verification

```
python3 verify/verify.py skew/1372
gzip -dc skew/1372/H1372.csv.gz | ./check_matrix
```

SHA-256 of `H1372.csv`:
`32401cf85e1d4b1284c1def607db2fa1e89e90d82babde46075c73478273f5c5`

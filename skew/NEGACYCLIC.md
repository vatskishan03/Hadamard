# Skew-Hadamard matrices of order 8v from negacyclic Goethals–Seidel arrays

Four orders, 1528 = 8 · 191, 2216 = 8 · 277, 2392 = 8 · 299 and 4712 = 8 · 589, come from a Goethals–Seidel array whose
blocks are negacyclic matrices of order 2v (v odd), searched in an equivalent quaternary form over ℤ_v. They improve
Cati–Pasechnik's Table 5 from 191(226), 277(5), 299(4) and 589(5) (589(4) through [9424](9424)) to 191(3), 277(3),
299(3) and 589(3).

| Order | v | Symmetry H | Orbits | Sums of the ±1 sequences | Before | Status |
|------:|---|---|---:|---|---|---|
| [1528](1528) | 191 | ⟨180⟩, order 19 | 11 | 6, −6, −22, 18 | 191(226) | New |
| [2216](2216) | 277 | ⟨273⟩, order 23 | 13 | −28, −34, −14, 32 | 277(5) | New |
| [2392](2392) | 299 = 13 · 23 | ⟨185, 209⟩, order 33 (Paley plug-in over ℤ₁₃ × F₂₃) | 15 | −60, −8, 12, 68 | 299(4) | New |
| [4712](4712) | 589 = 19 · 31 | ⟨311, 381⟩, order 45 (Paley plug-in over ℤ₁₉ × F₃₁) | 21 | −10, 34, −46, 12 | 589(5) | New |

## Construction

**The array.** Let n = 2v. For a ±1 sequence x of length n let N(x) be the negacyclic matrix with N[i, j] = x_{j−i}
for j ≥ i and −x_{n+j−i} for j < i, and let R be the back-diagonal identity. Negacyclic matrices commute, and
R N(x) R = N(x)ᵀ, so the Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

with A, B, C, D = N(a), N(b), N(c), N(d) is a Hadamard matrix of order 4n exactly when the negaperiodic
autocorrelations Σ_i x_i x_{i+s} (with x_{i+n} = −x_i) of the four sequences sum to 4n at s = 0 and to 0 elsewhere.
If a₀ = 1 and a_j = a_{n−j} for 0 < j < n, then A + Aᵀ = 2I and H is skew-Hadamard. Finlayson and Seberry
(Orthogonal designs from negacyclic matrices, Australas. J. Combin. 30 (2004) 319–330, Theorem 4) gave this array with
negacyclic blocks of odd order 3, 5 and 7, that is for orders 12, 20 and 28.

**The quaternary form.** Since v is odd, x ↦ i·g is an isomorphism ℤ[x]/(x^{2v} + 1) → ℤ[i][ℤ_v] (g a generator of
ℤ_v). A ±1 sequence of length 2v becomes (1 + i)·u with u a sequence on ℤ_v with values in {1, i, −1, −i}, the
negaperiodic condition becomes Σ u_X u_X* = 4v δ (complex periodic autocorrelations over ℤ_v), and the skew condition
on A becomes u(−g) = i · conj(u(g)) for g ≠ 0, with u(0) ∈ {1, −i}.

**Symmetry.** The search takes the four quaternary sequences constant on the orbits of a multiplier group H ≤ ℤ_v* of
odd order (so −1 ∉ H), with n_H orbits. A skew block then has 2^{n_H} choices and a free block 4^{n_H}, against
2^{(n_H − 1)/2} and 2^{n_H} for real skew Goethals–Seidel families with the same H; the number of equations doubles
(real and imaginary parts).

**Paley plug-ins.** For v = |K| q with q ≡ 3 (mod 4) a prime and gcd(|K|, q) = 1, take H = H_K × S, with H_K ≤ ℤ_|K|*
of odd order and S the squares of F_q*. A block is then s ⊗ δ + t ⊗ Q + e ⊗ J over K × F_q with Gaussian-integer
coefficients (x₀ at (k, 0), x₁ = e + t on the squares, x₂ = e − t on the non-squares, s = x₀ − e), exactly as in the
real plug-in of [search/NOTES.md](../search/NOTES.md) §15 and §87, and the compressions w = x₀ + (q − 1)e to K satisfy
Σ_X w_X w_X* = 4|K| q δ. At a position with x₂ = −x₁ (no defect) t is a free unit; a full defect (x₂ = x₁) and a half
defect (x₂ = ±i x₁) cost (q − 1) + 2 Re(x₀ ē) and (q − 1)/2 + 2 Re(x₀ ē), and the costs of all defects add up to
4|K|. The skew condition puts a half defect in A at k = 0. In 2392 and 4712 every defect is at k = 0: A's half defect
and two full defects (C and D, weights 10 + 20 + 22 = 52; B and D, 16 + 32 + 28 = 76).

## Search

`search/negaskew`: `enum_shapes.py` lists the multiplier groups of odd order of every target v with at most 13 orbits,
`qskew.c` joins the four orbit-level lists exactly (sums bucketed, hash join on the exact autocorrelation vectors),
`qp.c` is the two-level engine of the Paley plug-ins, and `certify.py` and `plugcert.py` turn a solution into the ±1
sequences and the certificate.

- *Targets.* The census-undecided v ≤ 2499 that the repository has not settled, for v ≤ 999 only those with a Table 5
  entry v(m), m ≥ 4.
- *Multiplier shapes.* 373 shapes with at most 13 orbits when the sweep ran (371 now, without 473 and 933): 2 hits
  (191 with 11 orbits, found in 4 s; 277 with 13 orbits, about 3 to 4 minutes), 370 EXHAUSTIVELY empty, and one
  (ℤ₂₁₉₇, 13 orbits) unfinished. The real class with the same H is empty at 191 and 277, as are all 510 real shapes
  with up to 17 orbits (508 now).
- *Plug-in shapes* (Table 5 targets, H_K with 5 to 9 orbits on K, so 15 to 27 orbits per block). Hits at 299 (ℤ₁₃,
  q = 23, H_K = ⟨3⟩: 13,056 solutions, exhaustive, 11 s) and 589 (ℤ₁₉, q = 31, ⟨7⟩: 3.3% of the level-1 instances
  complete). EXHAUSTIVELY empty: 539 (ℤ₄₉, q = 11, ⟨2⟩), 667, 707, 851, 899, 589 (ℤ₃₁, q = 19), 817 (two shapes),
  989 and 413. The defect costs cannot add up to 4|K| at 295, 335, 515, 611, 655, 695, 927 and 995. Sampled without a
  hit: 721, 737, 869, 973 and four shapes with 9 orbits on K; 597 and the other shapes with 9 orbits on K were not
  run.
- *Validation.* In real mode `qskew` and `qp` reproduce exact counts of real skew quadruples (176,256 for ℤ₃ × F₇,
  110,592 for ℤ₃ × F₁₁, 74,880,000 for ℤ₅ × F₇); `qp` and `qskew` give the same quaternary counts (79,164,864 for
  ℤ₃ × F₇, 3,753,216 for ℤ₃ × F₁₁, 3840 for ℤ₇ with ⟨2⟩ × F₁₁, normalised). Quaternary solutions over ℤ₃ and ℤ₅ give
  skew-Hadamard matrices of orders 24 and 40.

**One level up.** For v odd, ℤ[x]/(x^{4v} + 1) ≅ ℤ[ζ₈][ℤ_v], so negacyclic blocks of order 4v give skew orders 16v
from sequences on ℤ_v over a 16-letter alphabet. Its solutions do not all come from quadruples of the order-8v family
(at ℤ₆₁ with |H| = 15 it has solutions and the order-8v family has none), but a row-sum argument shows that a block
constant on the orbits of H with |H| = h needs (2h − 1)² + 15 ≤ 16v. All its multiplier shapes with at most 14 orbits
are therefore infeasible at the Table 5 entries with m ≥ 5
([search/NOTES.md](../search/NOTES.md) §89).

## Status

**New.** Table 5 gives 191(226) (from Paley's construction with a prime 2^{226} · 191 − 1), 277(5), 299(4) and 589(5).
The four v are also undecided in the census (order 4v), which these matrices do not settle. We found the negacyclic
Goethals–Seidel array for skew matrices only in Finlayson and Seberry (2004), with blocks of odd order, and the
quaternary multiplier search and the plug-ins nowhere. Doubling gives the skew orders 3056, 4432, 4784, 9424, …

## Verification

```
python3 verify/verify.py skew/1528 skew/2216 skew/2392 skew/4712
gzip -dc skew/1528/H1528.csv.gz | ./check_matrix
```

`verify/verify.py` (construction `goethals-seidel-negacyclic`) checks the negaperiodic condition, builds the array and
checks H Hᵀ = nI, H + Hᵀ = 2I and the SHA-256 of the canonical CSV.

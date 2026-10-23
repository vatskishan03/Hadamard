# Skew-Hadamard matrices of order 8v from negacyclic Goethals–Seidel arrays

Two orders, 1528 = 8 · 191 and 2216 = 8 · 277, come from a Goethals–Seidel array whose blocks are negacyclic matrices
of order 2v (v odd), searched in an equivalent quaternary form over ℤ_v. They improve Cati–Pasechnik's Table 5 from
191(226) and 277(5) to 191(3) and 277(3).

| Order | v | Symmetry H | Orbits | Sums of the ±1 sequences | Before | Status |
|------:|---:|---|---:|---|---|---|
| [1528](1528) | 191 | ⟨180⟩, order 19 | 11 | 6, −6, −22, 18 | 191(226) | New |
| [2216](2216) | 277 | ⟨273⟩, order 23 | 13 | −28, −34, −14, 32 | 277(5) | New |

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
(Australas. J. Combin. 30 (2004) 319–330) used negacyclic Goethals–Seidel arrays for orders 12, 20 and 28.

**The quaternary form.** Since v is odd, x ↦ i·g is an isomorphism ℤ[x]/(x^{2v} + 1) → ℤ[i][ℤ_v] (g a generator of
ℤ_v). A ±1 sequence of length 2v becomes (1 + i)·u with u a sequence on ℤ_v with values in {1, i, −1, −i}, the
negaperiodic condition becomes Σ u_X u_X* = 4v δ (complex periodic autocorrelations over ℤ_v), and the skew condition
on A becomes u(−g) = i · conj(u(g)) for g ≠ 0, with u(0) ∈ {1, −i}.

**Symmetry.** The search takes the four quaternary sequences constant on the orbits of a multiplier group H ≤ ℤ_v* of
odd order (so −1 ∉ H), with n_H orbits. A skew block then has 2^{n_H} choices and a free block 4^{n_H}, against
2^{(n_H − 1)/2} and 2^{n_H} for real skew Goethals–Seidel families with the same H; the number of equations doubles
(real and imaginary parts).

## Search

`search/negaskew`: `enum_shapes.py` lists the multiplier groups of odd order of every target v with at most 13 orbits,
`qskew.c` joins the four orbit-level lists exactly (sums bucketed, hash join on the exact autocorrelation vectors),
and `certify.py` turns a solution into the ±1 sequences and the certificate.

- *Targets.* The v ≤ 999 with Table 5 entries v(m), m ≥ 4, and the v ≤ 2499 undecided in the Đoković–Seberry census.
- *Result.* 373 shapes with at most 13 orbits: 2 hits (191 with 11 orbits, found in 4 s; 277 with 13 orbits, 3
  minutes), 370 EXHAUSTIVELY empty, and one (ℤ₂₁₉₇, 13 orbits) unfinished. The real class with the same H is empty at
  191 and 277, as are all 510 real shapes with up to 17 orbits.
- *Validation.* In real mode the same program reproduces exact counts of real skew quadruples (176,256 for ℤ₃ × F₇,
  110,592 for ℤ₃ × F₁₁, 74,880,000 for ℤ₅ × F₇); quaternary solutions over ℤ₃ and ℤ₅ give skew-Hadamard matrices of
  orders 24 and 40.

## Status

**New.** Table 5 gives 191(226) (from Paley's construction with a prime 2^{226} · 191 − 1) and 277(5). Both v are also
undecided in the census (order 4v), which these matrices do not settle. We found the negacyclic Goethals–Seidel array
for skew matrices only in Finlayson and Seberry (2004), at orders 12, 20 and 28, and the quaternary multiplier search
nowhere. Doubling gives the skew orders 3056, 4432, …

## Verification

```
python3 verify/verify.py skew/1528 skew/2216
gzip -dc skew/1528/H1528.csv.gz | ./check_matrix
```

`verify/verify.py` (construction `goethals-seidel-negacyclic`) checks the negaperiodic condition, builds the array and
checks H Hᵀ = nI, H + Hᵀ = 2I and the SHA-256 of the canonical CSV.

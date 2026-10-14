# Skew-Hadamard matrix of order 1900

1900 = 4 · 475 = 5 · 20 · 19, with 475 = 5² · 19.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I. This one comes from Seberry's theorem (Ars Combin. 6
(1978), Theorem 7): an orthogonal design of type (1, m, mn − m − 1) in order mn and amicable Hadamard matrices of order
n give a skew-Hadamard matrix of order mn(n − 1). Here n = 20 (a classical amicable order) and m = 5. The new
ingredient is an orthogonal design of type (1, 5, 94) in order 100.

**1. Amicable cores of order 19.** Paley's skew-Hadamard matrix S = I + C of order 20 (q = 19) and P = (−1) ⊕ (x ↦ −x)
satisfy PCP = −C, so xP + yCP is a symmetric orthogonal design of type (1, 19) and (S, P + CP) an amicable pair (the
classical amicable order q + 1, q ≡ 3 (mod 4)). Normalised at the fixed point of P with sign +1 (as in
[skew/2652](../2652)), they give W (order 19) skew with W1 = 0 and (I + W)(I + W)ᵀ = 20 I − J, and D symmetric with
D1 = −1 and D² = 20 I − J, with WD + DW = 0.

**2. A skew-type OD(100; 1, 5, 94).** Four circulant blocks of length t = 25 in the Goethals–Seidel array (codes in
`certificate.json`: 1 = a, ±2 = ±b, ±3 = ±c):

    A₁ = aI + cW₁,   A₂ = b(δ₀ + δ₁₀) + cW₂,   A₃ = b(δ₀ − δ₁₀) + cW₃,   A₄ = bδ₀ + cW₄,

with W₁ skew and Wᵢ = ±1 off the a and b positions. The b parts have autocorrelation sum
(2δ + δ±₁₀) + (2δ − δ±₁₀) + δ = 5δ. The array is an orthogonal design of type (1, 5, 94) exactly when

- Σᵢ PAF(Wᵢ)(s) = 0 for s ≠ 0 (the c² terms), and
- X(s) = W₂[s] + W₂[−s] + W₂[10 + s] + W₂[10 − s] + W₃[s] + W₃[−s] − W₃[10 + s] − W₃[10 − s] + W₄[s] + W₄[−s] = 0 for
  every s (the bc terms; the ac terms vanish because W₁ is skew).

Summing X over s gives 2r₂ + r₄ = 0 for the row sums. Together with 5r₂² + r₃² = 94 this forces
(r₂, r₃, r₄) = (3, 7, −6) up to signs.

**3. The plug-in.**

    H = I₁₀₀ ⊗ (I + W) + B ⊗ J + C ⊗ D        (order 100 · 19).

H + Hᵀ = 2I, the cross terms vanish as in [skew/7956](../7956), and

    H Hᵀ = I ⊗ [(I + W)(I + W)ᵀ + 5 J² + 94 D²] = I ⊗ [20 I − J + 5 · 19 J + 94 (20 I − J)] = 1900 I.

**Why m = 5 works in this form.** In the Goethals–Seidel form the b part of A₁ must be skew, so its row sum is 0, and
m = r₂(b)² + r₃(b)² + r₄(b)² must be a sum of three squares. 5 = 2² + 1² + 0² is. m = 7, 15, 23, … (≡ 7 mod 8) are
not, so the census values 1323, 759, 1815 and 2185 cannot be reached this way (`attack/NOTES.md` §31).

## Search

`attack/odm_sa.c` anneals the c parts with the b part fixed. The row sums of W₂, W₃, W₄ are held at their targets by
swapping a +1 and a −1 within a block, and pairs of the skew W₁ are flipped. The energy is
Σ (Σᵢ PAF(Wᵢ)(s))² + 2 Σ X(s)². The command is
`odm_sa 25 1 1000 1000000 40 1 2 "2:0:1,2:10:1,3:0:1,3:10:-1,4:0:1" 3 7 -6`. It found the design at restart 21, about
20 seconds on one core. The b templates with d = 5 and d = 1 in place of 10 gave nothing in 60 seconds.

## Reproduction

| | |
|---|---|
| Code | `attack/odm_sa.c` (search); `attack/sym_mm1.py` (`paley`, `skew_od_from_codes`, `check_skew_od`, `skew_seberry`); `verify/verify.py` (`build_skew_amicable_plug`, "paley" and "od" forms) |
| Command | `gcc -O2 -o odm_sa attack/odm_sa.c -lm && ./odm_sa 25 1 1000 1000000 40 1 2 "2:0:1,2:10:1,3:0:1,3:10:-1,4:0:1" 3 7 -6` |
| Randomness and caps | seed 1 (xorshift), 1000 restarts of 10⁶ steps, T from 40 to 1 |

`verify/verify.py` (construction `skew-amicable-plug`, "paley": 19) rebuilds Paley's design of order 20 and checks it.
It then checks the amicable pair and every identity of its cores, rebuilds the Goethals–Seidel design from the codes
and checks that it is a skew-type OD(100; 1, 5, 94) with 100 = 5 · 20. Finally it assembles H and checks
H Hᵀ = 1900 I, H + Hᵀ = 2I and the SHA-256 of the canonical CSV.

## Status

To our knowledge this is the first skew-Hadamard matrix of order 1900. The literature audit is in progress. Until it
is recorded here, read this as "found by us", not as "new". What has been checked so far:

1. **The Đoković–Seberry census** (PDF of 2025-05-28): 475 is in its list of 737 undecided skew values
   ("… 469, 473, 475, 479, …") and in its list of undecided symmetric values.
2. **Cati–Pasechnik** (arXiv:2411.18897), Table 5: "475(4)". The smallest known skew order 2ᵐ · 475 there is
   16 · 475 = 7600.
3. **Seberry–Yamada** (1992), Table 7.1: 475 has t = 4.
4. **SageMath** (develop branch of 2026-10-11, its existence logic traced, not executed): `Unknown`. Its
   orthogonal-design branch tries the split n = 20, m = 5 and needs an OD(100; 1, 5, 94), which is not in its table.
5. **Classical skew families** (arithmetic):
   - Paley needs 1899 = 3² · 211 to be a prime power.
   - 2(q + 1) needs q = 949 = 13 · 73.
   - 4(q + 1) needs q = 474.
   - Spence's q² + q + 1 skips 475 (463, 507).
   - Seberry's q(q + 1) and 2ʳn(2ʳn − 1) have no integer solution (43 · 44 = 1892).
   - (p − 1)^u + 1 needs 1899 to be a perfect power.
   - (|t| + 1)(q + 1) has no factorisation with q ≡ 5 (mod 8) a prime power.
   - Kronecker products and doubling are impossible (950 is not a Hadamard order).
   - The only splits of Seberry's theorem are n = 20, m = 5 (this one) and the degenerate n = 2.

The question is therefore whether an OD(100; 1, 5, 94) was known. We found none, and the census of 2025, by an author
of Theorem 7, lists 475 as undecided.

## Files

| File | Contents |
|---|---|
| `certificate.json` | q = 19 for the amicable pair ("paley") and the OD(100; 1, 5, 94) ("od": type and the four blocks as codes) |
| `H1900.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py skew/1900
gzip -dc skew/1900/H1900.csv.gz | ./check_matrix
```

SHA-256 of `H1900.csv`: `b607a5a5fb1a1a59024276b41c119684dca6f5807849e7ac5cfb34232dc531e2`

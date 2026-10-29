# Hadamard matrix of order 8696

8696 = 8 · 1087, with 1087 a prime ≡ 3 (mod 4).

## Construction

Seberry and Yamada's Hadamard matrices of order 8q (1992 survey, Lemma 8.12 and Corollaries 8.14 and 8.18), with
symmetric Williamson-type matrices of order n = (q − 1)/2 = 543.

**1. Williamson-type matrices of order 543.** 541 is a prime ≡ 5 (mod 8) and (541 + 1)/2 = 271 is a prime. So Seberry
and Yamada's Corollary 27(i) (JCMCC 7, 1990) gives four symmetric ±1 matrices W₁, …, W₄ of order 541 + 2 = 543 that
commute pairwise, with Σ W_k² = 4 · 543 I. It is their Lemma 26 with U = (I, 0, BR, BR) and V = (M, N, R, R), where B is
the Paley core of order 271, R the reversal, and M and N are Yamada's symmetric circulants of order 271 with
M² + N² = 541 I (`search/symwt.py "cor27i(541)"`). The certificate stores the four matrices, and `verify/verify.py`
checks the conditions exactly.

**2. The 8q construction.** Let q = 1087, g = 3 a primitive root, s_a = g^(2a), χ the Legendre symbol,
A = [χ(s_a − s_b)] (skew, row sums 0), B = [χ(s_a + s_b)] (symmetric, row sums −1), so that B² − A² = qI − 2J, and R
the permutation a ↦ −a (mod n). Put

```
U = (−BR, −BR, A, A, 0, 0, 0, 0),   V = (0, 0, I, I, W₁, W₂, W₃, W₄),
X_j = [[1, s_j eᵀ], [s_j e, U_j ⊗ J₂ + V_j ⊗ [[1, −1], [−1, 1]]]],   s_j = −1 (j = 1, 2), +1 (j = 3, …, 8).
```

These are eight ±1 matrices of order q. Because J₂ [[1, −1], [−1, 1]] = 0, the products split into a Paley part and a
Williamson-type part, and the X_j are pairwise amicable with Σ X_j X_jᵀ = 8q I. Plugged into an orthogonal design
OD(8; 1⁸), they give the Hadamard matrix of order 8 · 1087 = 8696. `verify/verify.py` checks the amicability and the sum
exactly before plugging, and then H Hᵀ = 8696 I.

## Search

None. The route came from a closure computation over the objects of this repository under the recursive rules of
the tables ([search/NOTES.md](../../search/NOTES.md) §95). The 8q theorem had been set aside for q = 1087 in the
analysis of the orders 8p (§92), because the 2026 survey excludes Williamson-type matrices of order 543.
`search/symwt.py` already built them.

**Reproduction.**

| | |
|---|---|
| Code | `search/sy8q.py`, `search/symwt.py`, `verify/verify.py` (construction `seberry-yamada-8q`) |
| Command | `cd search; python3 sy8q.py wt 1087 "cor27i(541)"` |
| Randomness and caps | none, deterministic |
| Result | `certificate.json`; matrix SHA-256 `81eafe79…` |

Two independent rebuilds with their own code checked the matrix. One reproduced it byte for byte, the other used another
primitive root and a different matrix; both verified H Hᵀ = 8696 I exactly.

## Status

**Known.** The order was recorded as known through a valid route that the 2026 survey dropped. The Hadamard Matrix
Explorer (math74.com, version 53, October 2026, not dated more precisely) also posts a matrix of this order built by the
same route.

- Seberry–Yamada's Table A.1 (1992) gives Williamson-type matrices of order 543 with key wi: order q + 2 for a prime
  power q ≡ 1 (mod 4) with (q + 1)/2 a prime power (Miyamoto), here q = 541 and (q + 1)/2 = 271, which is their
  Corollary 27(i). With them their Table A.2, and Cati–Pasechnik's Table 4 (1087(3)), record order 8 · 1087.
- The survey of Egan, Kotsireas, London and Ó Catháin (arXiv:2610.27175, October 2026) takes Williamson-type matrices
  from Table A.1 "except 543 and those resting on its Corollary 8.8(1)" and lists 8696 among its six open orders 8p up
  to 10 000 (Table 5). The key-wi route does not rest on Corollary 8.8(i), and it works: the matrices of order 543 are
  built and checked here, and so is the matrix of order 8696. The survey's open list should therefore drop 8696.
- No other route that we know:
  - Paley's constructions need 8695 = 5 · 37 · 47 or 4347 = 3³ · 7 · 23 to be prime powers.
  - Doubling needs a Hadamard matrix of order 4348 = 4 · 1087, which is not known ([FRONTIER.md](../../FRONTIER.md)).
  - Miyamoto's theorem for order 8q needs 2q − 3 = 2171 = 13 · 167 to be a prime power.
  - The complex Golay decomposition 8(x + y) needs two complex Golay lengths with sum 1087.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the four Williamson-type matrices (packed upper triangles), the SHA-256 of the matrix |

The matrix (8696² entries) is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 8696
```

SHA-256 of `H8696.csv`: `81eafe79dcfb1b945fef22652f6f963748a84c0d8325249a8d4d42877c6fbc5f`

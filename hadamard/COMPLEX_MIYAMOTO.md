# Miyamoto's Theorem 7 in quaternary form

Miyamoto's theorem (J. Combin. Theory Ser. A 57 (1991) 86–108) builds a Hadamard matrix of order 4q from Paley's
symmetric conference matrix of order q + 1 and a Hadamard matrix of order q − 1, for a prime power q ≡ 1 (mod 4). For
q ≡ 3 (mod 4) Paley's conference matrix is skew and there is no Hadamard matrix of order q − 1 ≡ 2 (mod 4). Both
obstacles disappear over the units {1, −1, i, −i}: i times a skew conference matrix is Hermitian, and quaternary
Hadamard matrices of order 2N, N odd, exist for many N ([complex-hadamard](../complex-hadamard/README.md)).

**Theorem.** Let q ≡ 3 (mod 4) be a prime and K a quaternary Hadamard matrix of order q − 1. Then there is a quaternary
Hadamard matrix of order 4q, and by Turyn's doubling a Hadamard matrix of order 8q.

**Credit.** This is very probably Miyamoto's own Theorem 7, stated here with the intermediate quaternary matrix:
- His §6, "Hadamard matrices of order 8q and 16q", uses complex Hadamard matrices to obtain orders 8q.
- The Hadamard Matrix Explorer (version 53, notation page) describes Theorem 7 as "a complex Hadamard matrix of order
  q − 1 … gives order 8q".
- Seberry and Yamada (1989, 1992, key m2), Brent–Osborn–Smith (2013) and Egan, Kotsireas, London and Ó Catháin (2026,
  "Miyamoto III") cite it only when 2q − 3 is a prime power, where Paley's matrix of order 2q − 2 supplies the complex
  input.

We could not read §6 ([NOVELTY.md](../NOVELTY.md)). The construction and proof below are our own and self-contained.
`verify/verify.py` implements them (`complex_miyamoto`, certificate type `complex-miyamoto`). They are checked exactly:
- for q = 7, 19, 31 and 43 with K from Paley's complex construction;
- for q = 223, 263, 331 and 431 with the matrices of orders 222, 262, 330 and 430 of
  [complex-hadamard](../complex-hadamard/README.md);
- for q = 1559 ([12472](12472)).

## Construction

Let m = (q − 1)/2, χ the quadratic character of GF(q), and order GF(q)* as the m squares followed by the m non-squares.

- *The conference part.* W = i [χ(y − x)] on GF(q)* (order 2m) is Hermitian with zero diagonal, since χ(−1) = −1.
  (It is the core of i times Paley's skew conference matrix on ∞, 0 and GF(q)*.)
- *The Hadamard part.* V = [[I, K], [−K*, I]] (order 4m).
- U = diag(W, W) (order 4m). U and V are cut into 4 × 4 blocks U_rc, V_rc of order m. Their supports are complementary:
  V is nonzero exactly on the diagonals of the diagonal blocks and on the off-diagonal blocks of order 2m.
- *Borders.* B is the Sylvester matrix of order 4, F = [[1, 1, 1, 1], [1, 1, −1, −1], [1, −1, 1, −1], [−1, 1, 1, −1]],
  P the permutation (0 1)(2 3), R = [[0, i, 0, 0], [−i, 0, 0, 0], [0, 0, 0, i], [0, 0, −i, 0]],
  B′ = R P F and A = −B P F/2. A, B and B′ are quaternary Hadamard matrices of order 4.

H is the 4 × 4 array of blocks of order q = 1 + 2m

```
H_rc = [ A_rc        B_rc e                               ]
       [ B′_rc eᵀ    [[U_rc + V_rc, U_rc − V_rc],
                      [U_rc − V_rc, U_rc + V_rc]]           ],
```

with e the all-ones row of length 2m. Its entries are 1, −1, i, −i.

## Proof

Index the rows of block row r by ∞ and by (h, x), h ∈ {0, 1}, x ∈ [m]. Three facts about the ingredients:

1. V V* = [[I + K K*, 0], [0, K* K + I]] = qI.
2. W W* = qI − J − χχᵀ, χ = (1, …, 1, −1, …, −1). This holds because Σ_{z ∈ GF(q)} χ(z − x)χ(z − y) = qδ − 1 and the
   term z = 0 contributes χ(x)χ(y). So W W* has blocks qI − 2J on the diagonal and 0 off it.
3. The rows of each block of W have constant sums. For x a square, Σ over squares y of χ(y − x) is 0 and over
   non-squares is 1. For x a non-square these are −1 and 0. Both follow from Σ_{y ≠ 0} χ(y − x) = χ(x) and
   Σ_{y ≠ 0} χ(y − x)χ(y) = −1. So the block row sums of U form the matrix R above (times i), a monomial matrix with
   unit entries.

In the inner product of two core rows (r, h, x) and (r′, h′, x′), the U and V terms do not mix, because
J₂ ⊗ U + D₂ ⊗ V with D₂ = [[1, −1], [−1, 1]] has J₂D₂ = 0. It equals

    Σ_c B′_rc conj(B′_r′c) + 2 (Σ_c U_rc U_r′c*)[x, x′] + 2(−1)^(h+h′) (Σ_c V_rc V_r′c*)[x, x′].

By fact 1 the V term is 2q δ. By fact 2 the U term is 2(qI − 2J) on r = r′ and 0 otherwise, so the sum is 4q δ exactly
when B′B′* = 4I.

The border row (r, ∞) against a core row (r′, h′, x′) gives Σ_c A_rc conj(B′_r′c) + 2 Σ_c B_rc conj(R_r′c), since the V
terms cancel over h. So we need A B′* + 2 B R* = 0. Two border rows give A A* + 2m B B* = 4q I. With B′ = R P F and
A = −B P F/2:
- R* R = P P* = I and F F* = 4I, so B′B′* = R P F F* P* R* = 4I and A A* = B P F F* P* B*/4 = B B* = 4I;
- A B′* = −B P F F* P* R*/2 = −2 B R*.

So H H* = 4qI. The construction is the M-structure of Miyamoto's theorem (Seberry–Yamada's layout in `verify.py`'s
`miyamoto`) with a Hermitian core in place of the symmetric one. The difference is fact 3: the row sums of the
Hermitian core sit on the opposite halves, so the border is R P F instead of a sign pattern.

## Consequences

Turyn's doubling of the quaternary matrix of order 4q gives a Hadamard matrix of order 8q. So every prime q ≡ 3 (mod 4)
with a known quaternary Hadamard matrix of order q − 1 has a Hadamard matrix of order 8q. The tables use the theorem
only when 2q − 3 is a prime power; the quaternary matrices of 2026 widen it. Against Table 4 of Cati and Pasechnik
(arXiv:2411.18897v2; entries n(m) with m ≥ 4, n prime ≡ 3 mod 4):

| q | Table 4 | quaternary matrix of order q − 1 | result |
|---:|---|---|---|
| 1559 | 1559(4) | 1558, constructed here and on Zaitsev's page ([C1558.json](../complex-hadamard/C1558.json)) | order 8 · 1559 = 12472 ([12472](12472)); the quaternary matrix of order 6236 = 4 · 1559, where a real Hadamard matrix is open |
| 491 | 491(5) | 490, on Zaitsev's page | order 8 · 491 = 3928, already known by doubling: 491 is no longer open ([FRONTIER.md](../FRONTIER.md)) |

For the other such q the quaternary matrix of order q − 1 is unknown. That includes the three open orders 8p of Egan,
Kotsireas, London and Ó Catháin: 4792, 5752 and 5912 would follow from quaternary Hadamard matrices of orders 598,
718 and 738, all open ([complex-hadamard/README.md](../complex-hadamard/README.md)).

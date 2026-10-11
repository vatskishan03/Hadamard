# Method

The Goethals–Seidel orders in this repository came from a structured search for
supplementary difference sets in ℤ_p. This note describes the search space, the
algorithm, and the cost model used to choose targets.

## Search space

Let p be prime and let H ≤ ℤ_p* be a subgroup of odd order h. A ±1 sequence on
ℤ_p that is constant on the cosets of H is fixed by one sign per coset together
with its value at 0, so it has (p − 1)/h + 1 free signs. Its periodic
autocorrelation P(s) = Σₜ x(t) x(t + s) is constant on the classes ±sH. There are
C = (p − 1)/(2h) such classes.

Four such sequences give a Hadamard matrix of order 4p through the
Goethals–Seidel array exactly when

```
P_A(s) + P_B(s) + P_C(s) + P_D(s) = 0                                         (1)
```

for one representative s of every class. This is a system of C integer equations,
of which C − 1 are independent: for each sequence with row sum r, the class values
weighted by class size add up to r² − p.

The row sums satisfy r ≡ ±1 (mod 2h). This restricts the admissible patterns
r_A² + r_B² + r_C² + r_D² = 4p, and for many pairs (p, h) no pattern exists.

## Algorithm

**Candidate generation.** For each required row sum, a TPU generates 10⁸ to 10⁹
sequences by batched local search. The search minimises Σ_c (P(c) − m_c)², where
m_c is the exact mean of P(c) over all sequences with that row sum. The moves swap
the signs of two cosets, which preserves the row sum, and iterated local search
adds perturbation and re-descent. This concentrates the class values near their
means. For p = 523 and h = 9 the per-class variance falls from about 490 to 66.

**Matching.** With v_c = (P(c) − (p mod 4))/4, condition (1) reads

```
v_A + v_B + v_C + v_D = −(p mod 4)
```

on every class. We use the extended k-tree (generalised birthday) algorithm of
Wagner [1] and Minder–Sinclair [2]:

1. Fix the partial sums τ of v_A + v_D on s₁ classes.
2. Enumerate the pairs (A, D) and (B, C) whose partial sums match, by merging
   lists sorted on additive 64-bit hashes.
3. Join the two pair lists on a hash of the remaining classes.

The values τ are processed in decreasing order of probability, and every match is
checked exactly. With lists of 10⁸ elements, one value of τ takes a few seconds
on 44 cores.

**Multiplier images.** If x is constant on the cosets of H, then so is
x_u(t) = x(ut) for any unit u, and its class values are a permutation of those of
x. Every candidate therefore yields up to C candidates of the same quality at no
generation cost (`search/expand.c`). In practice the deepest candidates for a
given structure occupy a small region that is nearly closed under this action.
For p = 919, five images per candidate gave 1.8 times as many distinct sequences,
and the expected yield per τ rose by a factor of 39.

## Cost model

Fix a structure (p, h and a row-sum pattern). Let N be the number of sequences
with the required row sum and E the expected number of solutions. We estimate E
from the empirical class-value distributions, taking into account the dependent
class and the lattice of attainable values. Lists of size L selected by descent
then contain on the order of

```
E · (L/N)³ · (6/π)^((C−1)/2)
```

solutions. The work per solution found therefore grows like N³/(E·L): a
structure is easy when its pools are small and its solutions numerous. The model
matched the observed effort for order 2092, and we used it to rank every
unresolved entry of the Cati–Pasechnik table.

One family stands out. For primes p = m² + m + 1, a subgroup of order m or m + 1
makes the row sums (1, 1, 1, 2m + 1) admissible automatically. The pools are about
2³⁵ sequences, there are about 10¹² solutions, and the orders
[4492](4492) (m = 33) and [5932](5932) (m = 38) were found within seconds of
matching. Spence's construction covers the same orders only when m is a prime
power.

## Limitations

The method needs a multiplier subgroup with admissible row sums and a moderate
number of cosets. For prime p this is a real constraint.

**The √p bound.** Every row sum satisfies r ≡ ±1 (mod 2h). The pattern
(±1, ±1, ±1, ±1) is excluded because 4p > 4, so some |r| is at least 2h − 1. Then
(2h − 1)² ≤ 4p, which forces h ≤ √p + 1/2 and at least √p − 1 cosets. Pools
therefore grow like 2^√p, while the class values spread more as p grows. Above
p ≈ 2000 no structure is cheap. For [4412](4412), a pilot measured 10⁻⁹ expected
solutions per τ.

For [2148](2148) no usable structure exists at all.

**Small structures.** Structures too small for the Gaussian model to be trusted
were checked exhaustively instead. There were 1,319 such structures across the
unresolved entries; 1,159 were searched completely, and none contains a
solution. See [FRONTIER.md](FRONTIER.md).

We also tested generating pairs of sequences jointly so as to concentrate their
pair sums. It reduces the pair variance by a factor of about 2.5, but it cannot
compete with the combinatorial pairing of the k-tree.

## References

1. D. Wagner, A generalized birthday problem, CRYPTO 2002, LNCS 2442, 288–304.
2. L. Minder, A. Sinclair, The extended k-tree algorithm, J. Cryptology 25 (2012) 349–382.
3. D. Ž. Đoković, O. Golubitsky, I. S. Kotsireas, Some new orders of Hadamard and
   skew-Hadamard matrices, J. Combin. Des. 22 (2014) 270–277.

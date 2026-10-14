# Twisted Szöllősi quadruples

Two specialisations of Szöllősi's array (arXiv:2603.09572, Theorem 4). Each one turns the four-sequence search for a
complex Hadamard matrix of order 2N into a two-list search. To our knowledge neither has been stated before. They
are not new matrices by themselves; they are new families in which to look for them. Code: `attack/tsz.py` and
`attack/tsz.c`.

## Setting

Szöllősi's theorem needs four ±1 sequences A, B, C, D on ℤ_N such that:

- A and B are symmetric;
- Σ_k PAF_k(s) = 0 for every s ≠ 0, where PAF_X(s) = Σ_t X(t) X(t + s).

Its proof uses only that the group is abelian, so it holds verbatim over any abelian group G of order N.

Notation:

- H is a group of automorphisms of G (multipliers on ℤ_N) with −1 ∉ H, and L = ⟨H, −1⟩.
- μ is an automorphism with μ ∉ L and μ² ∈ L. L′ = ⟨L, μ⟩ contains L with index 2.
- X∘μ denotes the sequence t ↦ X(μt).

For every sequence X and automorphism ν, PAF_{X∘ν}(s) = PAF_X(νs), by the substitution t ↦ ν⁻¹t. If X is
constant on the orbits of a group K of automorphisms, PAF_X is constant on the K-orbits of G.

## Full twist

**Construction.** A is constant on the L-orbits, C is constant on the H-orbits, and B = A∘μ, D = C∘μ.

Then B is symmetric because A is: B(−t) = A(−μt) = A(μt). Put F = PAF_A + PAF_C. Then

    Σ PAF (s) = F(s) + F(μs).

F is constant on the L-orbits and μ² ∈ L, so this sum is constant on the L′-orbits. The condition therefore has to be
checked on one representative of each L′-orbit of G ∖ {0}, about half as many classes as L gives. Summing it over
all s ≠ 0 gives 2 r_A² + 2 r_C² = 4N, so

    r_A² + r_C² = 2N,

and N must be a sum of two squares. On ℤ_p with p prime this means p ≡ 1 (mod 4).

**Parity.** PAF_X(s) ≡ N (mod 4) for every ±1 sequence X of odd length N, so F ≡ 2 (mod 4). Suppose μ maps the
L-orbit of some s ≠ 0 to itself. Then the condition there reads 2F(s) = 0, which is impossible. Every such
structure is empty. This cannot happen on ℤ_p, but it can on ℤ_N with N composite.

## Half twist

**Construction.** A and B are constant on the L′-orbits (so both are symmetric), C is constant on the H-orbits, and
D = C∘μ. Then

    Σ PAF (s) = PAF_A(s) + PAF_B(s) + PAF_C(s) + PAF_C(μs),

which is again constant on the L′-orbits. The row sums satisfy r_A² + r_B² + 2 r_C² = 4N. No parity obstruction
arises: on an L-orbit that μ fixes, the last two terms are 2 PAF_C(s), and the sum is ≡ 4N ≡ 0 (mod 4).

[L′ : L] = 2 needs an element of order 2 in the unit group modulo L. When ℤ_N^* is cyclic of order 2 · odd, as for
p ≡ 3 (mod 4) prime, there is none, so the half twist is unavailable there. Every N with two distinct prime factors
has such twists.

**Bordered half twist.** The same holds for the bordered array of [BORDERED.md](BORDERED.md): sequences on a group
of order N − 1, Σ PAF = −4 off zero, row sums (r_A, r_B, r_C) = (0, 2, 0) up to swapping A and B. The full twist
cannot be bordered: it forces r_A = r_B, and the row sums of the bordered array are (±2, 0, 0, 0) up to order.

## Search

Both twists are two-list problems. On one side is the hashed list of A (full twist) or of the pairs (A, B) (half
twist). On the other side every C of the given row sum is enumerated by a depth-first search over the H-orbits, with
class values updated incrementally, and looked up. The cost is about the size of the C pool, so structures with up
to about 40 H-orbits are searched exhaustively in minutes. The Gaussian estimate is about 2^(1.5 O) · 2^(−b O/4), with
O the number of H-orbits and b ≈ 4 bits per class condition. That is the square root of the untwisted estimate for
the same H, at the cost of one pool instead of a four-list join.

## Verification

`python3 attack/tsz.py test` searches small cases exactly and rebuilds every hit as a complex Hadamard matrix
(K K* = 2N I in Gaussian integers) and as a Goethals–Seidel matrix of order 4N. Solution counts:

| N | kind | Gaussian estimate | solutions |
|---:|---|---:|---:|
| 13 | full, H = 1, μ = 5 | 10^3.15 | 780 |
| 17 | full, H = 1, μ = 4 | 10^3.87 | 6256 |
| 29 | full, H = 1, μ = 12 | 10^5.91 | 989016 |
| 41 | full, H = ⟨10⟩, μ = 2 | 10^1.76 | 48 |
| 61 | full, H = ⟨13⟩, μ = 11 | 10^3.57 | 4680 |
| 15 | half, H = 1, μ = 4 | 10^2.14 | 480 |
| 21 | half, H = 1, μ = 8 | 10^2.61 | 1008 |
| 39 | half, H = ⟨25⟩, μ = 5 | 10^1.01 | 144 |
| 51 | half, H = ⟨16⟩, μ = 4 | 10^1.41 | 288 |
| 27 | bordered half, H = 1, μ = 5 | 10^4.42 | 406224 |
| 31 | bordered half, H = 1, μ = 11 | 10^3.83 | 61440 |

For prime N the counts agree with the estimates within a factor of two. For composite N they exceed them, by up to
a factor of 15 here.

**A small family.** Let p ≡ 5 (mod 8) and let H be the fourth powers. Then L is the squares, μ is a non-square, and
L′ = ℤ_p^* is transitive on ℤ_p ∖ {0}, so Σ PAF is constant off zero. The full twist is then a solution exactly
when r_A² + r_C² = 2p. The orbit structure allows r_A ∈ {±1, ±(p − 2), ±p} and r_C ∈ {±1, ±1 ± (p − 1)/2,
±1 ± (p − 1)}. The equation then has solutions only for p = 5 (|r_A|, |r_C| = 1, 3 or 3, 1) and p = 13
(|r_A| = 1, |r_C| = 5), so this is not an infinite family. A direct check confirms that at p = 5 and p = 13 every
pair (A, C) with these row sums is a solution: 80 and 32 pairs.

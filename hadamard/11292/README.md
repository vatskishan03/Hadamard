# Hadamard matrix of order 11292

11292 = 4 · 2823, with 2823 = 3 · 941.

## Construction

Miyamoto's Hadamard matrix of order 3764 (q = 941, a prime ≡ 1 mod 4) is built from Paley's conference matrix of order
942 and a Hadamard matrix K of order 940. Here K is the Williamson array of four symmetric, pairwise commuting ±1
matrices of order w = 235 (symmetric Williamson-type matrices), which come from a symmetric conference matrix of order
118 found here. The four quadrants of K then lie in a commutative algebra, so Miyamoto's matrix is a 4 × 4 array of
bordered blocks of order 941 that commute and are transposed by one permutation. Multiplying it by T-matrices of order
3 in the Goethals–Seidel array gives the Hadamard matrix of order 4 · 941 · 3 = 11292. The construction and its proof
are in [MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#symmetric-williamson-type-matrices).

- Paley's conference matrix of order 942, its points listed as ∞, 0, g^(2j), g^(2j+1) (j < 470) with the primitive root
  g = 2.
- The symmetric conference matrix C of order 118 of [symmetric/940](../../symmetric/940)
  (`search/data/conference118_5fe566972f21790b.txt`, normalised: first row and column +), found by an exact cover
  ([search/NOTES.md](../../search/NOTES.md) §95). Its core Q of order 117 satisfies Q = Qᵀ, Q² = 117I − J.
- Symmetric Williamson-type matrices of order 235 = 233 + 2 by Seberry and Yamada's Corollary 27 (J. Combin. Math.
  Combin. Comput. 7 (1990)): their Lemma 26 with U = (I, 0, Q, Q) and V = (M, N, I, I), where M and N are Yamada's
  symmetric circulants of order 117 with M² + N² = 233I (over GF(233²)). Two of the four matrices are equal. K is their
  Williamson array [[X₁, X₂, X₃, X₄], [−X₂, X₁, −X₄, X₃], [−X₃, X₄, X₁, −X₂], [−X₄, −X₃, X₂, X₁]], relabelled so that
  the half swap of each quadrant becomes j ↦ 1 − j on ℤ₄₇₀. `search/symwt.py` builds them (recipe
  `cor27c(233, "conference118_5fe566972f21790b.txt")`) and checks that they are symmetric, commute and have
  X₁² + X₂² + X₃² + X₄² = 940I.
- T-matrices of order 3 (first rows, from base sequences BS(2, 1), built by the program): `+00`, `0+0`, `00+`, `000`.

## Search

The conference matrix of order 118 came from an exact cover ([search/NOTES.md](../../search/NOTES.md) §95); the rest
is determined by q, the recipe and the T-matrices.

| Code | Command | Randomness and caps | Result |
|---|---|---|---|
| `search/miyamoto_blocks.py`, `search/symwt.py`, `verify/verify.py` construction `miyamoto-blocks` | `python3 search/miyamoto_blocks.py certify 941 3 'wt:cor27c(233, "conference118_5fe566972f21790b.txt")'` | none | `certificate.json`, matrix SHA-256 `03b29b37…` |

## Status

**New.** Cati–Pasechnik's Table 4 lists n = 2823 as open, and we found no later construction.

- [MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#reach) wrote down this route before the conference matrix was known to us.
  It needs symmetric Williamson-type matrices of order 235 or a negaperiodic Golay pair of length 470. Corollary 27
  gives the matrices from a symmetric conference matrix of order (233 + 3)/2 = 118, whose existence was open until
  Brouwer's table of strongly regular graphs (read on 17 October 2026) cited J. Y. Z. Ho (email, October 2026) for the
  conference graph SRG(117, 58, 28, 29). Ho's matrix is not public. Ours was found independently and is Concurrent with
  it ([symmetric/940](../../symmetric/940)). We know of no construction of order 11292 by Ho or anyone else.
- Paley's constructions do not apply: 4n − 1 = 11291 = 7 · 1613 and 2n − 1 = 5645 = 5 · 1129 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4).
- Yamada–Kiyasu does not apply: n − 2 = 2821 = 7 · 13 · 31 is not a prime power.
- Cooper–Wallis needs Williamson-type matrices of order 941 or T-matrices of order 941. Before this construction
  neither was known to us. The repeated block of the matrices of order 235 now gives symmetric Williamson-type
  matrices of order 941 = 4 · 235 + 1 by Corollary 8.8(ii) (`symwt.py 'cor29r(941, cor27c(233, …))'`, checked), so
  Cooper–Wallis gives the order too, from the same conference matrix.
- The Paley-block route does not reach it: for the divisors s of 2823, 4s − 1 = 3763 = 53 · 71 and 11291 are not
  prime powers, and s = 1 and 3 need T-matrices of order 2823 or 941, which are not known to us.
- No negaperiodic or periodic Golay pair of length 470 is known to us; a periodic one cannot exist (940 = 2² · 5 · 47
  is not a sum of two squares).

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the four symmetric Williamson-type matrices (upper triangles, base64), the T-matrices and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it, after checking that the four matrices are symmetric ±1
matrices that commute pairwise and have X₁² + X₂² + X₃² + X₄² = 4wI.

## Verification

```
python3 verify/verify.py hadamard/11292
```

SHA-256 of `H11292.csv`: `03b29b373d81e48a7d7b0fbc180988342b8bb63f78e1d7c25d31a2510d1fea53`

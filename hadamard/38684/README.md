# Hadamard matrix of order 38684

38684 = 4 · 9671, with 9671 = 19 · 509.

## Construction

Miyamoto's Hadamard matrix of order 2036 (q = 509, a prime ≡ 1 mod 4) is built from Paley's conference matrix of order
510 and a Hadamard matrix K of order 508. Here K is the Williamson array of four pairwise commuting symmetric matrices
of order w = 127, obtained from near-Williamson matrices. Its four quadrants then lie in a commutative algebra, so
Miyamoto's matrix is a 4 × 4 array of bordered blocks of order 509 that commute and are transposed by one permutation.
Multiplying it by T-matrices of order 19 in the Goethals–Seidel array gives the Hadamard matrix of order 4 · 509 · 19 =
38684. The construction and its proof are in [MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#near-williamson-matrices).

- Paley's conference matrix of order 510, its points listed as ∞, 0, g^(2j), g^(2j+1) (j < 254) with the primitive root
  g = 2.
- Near-Williamson matrices of order 127 from the database of arXiv:2605.08661v2 (cs.uleth.ca/~zaitsevv/nw/NW127.txt),
  [near-williamson/NW127.txt](../../near-williamson/NW127.txt): circulant ±1 matrices A, B, C, D with B, C, D symmetric
  and AAᵀ + BBᵀ + CCᵀ + DDᵀ = 508I. With R the back-identity, AR, B, C, D are symmetric and commute pairwise, and K is
  their Williamson array [[X₁, X₂, X₃, X₄], [−X₂, X₁, −X₄, X₃], [−X₃, X₄, X₁, −X₂], [−X₄, −X₃, X₂, X₁]], relabelled so
  that the half swap of each quadrant becomes j ↦ 1 − j on ℤ_254.
- T-matrices of order 19 (first rows, [t-matrices/T19.txt](../../t-matrices/T19.txt), from `search/tmat_mult.py 19 7`):
  `++00000+-00+-00000-`, `0000000000000000000`, `00++0000000000+0000`, `0000+++00+-00-0-++0`.

## Search

None: the matrix is determined by q, the near-Williamson matrices and the T-matrices. `python3 search/miyamoto_blocks.py
certify 509 19 nw:near-williamson/NW127.txt t-matrices/T19.txt` writes the certificate.

## Status

**First explicit construction** ([MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#near-williamson-matrices)). Table 1.53 of
the Handbook of Combinatorial Designs (2nd ed., 2007) gives n = 9671 with t = 3. Đoković (arXiv:1008.2043, Propositions
3.2–3.4; Int. Math. Forum 6 (2011) 3061–3067) settles it with Williamson-type matrices of order 509, which he obtains
from his good matrices of order 127 through Seberry–Yamada's Corollary 8.8(i). The printed proof of that corollary does
not give Williamson-type matrices ([NOTES.md](../../NOTES.md) §7), and we found no other construction of the order. Good
matrices of order 127 also satisfy the hypotheses used here, so Đoković's matrices would serve as well as the
near-Williamson matrices of the 2026 database.

- Paley's constructions do not apply: 4n − 1 = 38683 = 101 · 383 and 2n − 1 = 19341 = 3² · 7 · 307 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4).
- Yamada–Kiyasu does not apply: n − 2 = 9669 = 3 · 11 · 293 is not a prime power ≡ 1 (mod 4).
- Cooper–Wallis needs Williamson-type matrices of order 509 or T-matrices of order 509; neither is known to us.
- The Paley-block route does not reach it: 4s − 1 = 75, 2035 and 38683 are not prime powers ≡ 3 (mod 8) for the other
  divisors s of 9671, and the remaining ones need T-matrices of order 9671 (s = 1), which are not known to us.
- The Miyamoto-block route with a negaperiodic or periodic Golay pair of length 254 does not apply: no negaperiodic pair
  of that length is known (Balonin–Đoković, list (31), quasi-Williamson matrices of order 127), and a periodic one
  cannot exist (508 = 2² · 127 is not a sum of two squares).
- SageMath 10.8 (passagemath 10.8.13, October 2026) does not construct it: `hadamard_matrix(38684, existence=True)`
  returns Unknown.

We found no valid construction of order 38684 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the near-Williamson matrices, the T-matrices and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 38684
```

SHA-256 of `H38684.csv`: `38b3fc09630077072b4af85c50643974ee311e3f6e72acf6072e8c3bcc59817c`

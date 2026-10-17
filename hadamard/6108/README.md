# Hadamard matrix of order 6108

6108 = 4 · 1527, with 1527 = 3 · 509.

## Construction

Miyamoto's Hadamard matrix of order 2036 (q = 509, a prime ≡ 1 mod 4) is built from Paley's conference matrix of order
510 and a Hadamard matrix K of order 508. Here K is the Williamson array of four pairwise commuting symmetric matrices
of order w = 127, obtained from near-Williamson matrices. Its four quadrants then lie in a commutative algebra, so
Miyamoto's matrix is a 4 × 4 array of bordered blocks of order 509 that commute and are transposed by one permutation.
Multiplying it by T-matrices of order 3 in the Goethals–Seidel array gives the Hadamard matrix of order 4 · 509 · 3 =
6108. The construction and its proof are in [MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#near-williamson-matrices).

- Paley's conference matrix of order 510, its points listed as ∞, 0, g^(2j), g^(2j+1) (j < 254) with the primitive root
  g = 2.
- Near-Williamson matrices of order 127 from the database of arXiv:2605.08661v2 (cs.uleth.ca/~zaitsevv/nw/NW127.txt),
  [near-williamson/NW127.txt](../../near-williamson/NW127.txt): circulant ±1 matrices A, B, C, D with B, C, D symmetric
  and AAᵀ + BBᵀ + CCᵀ + DDᵀ = 508I. With R the back-identity, AR, B, C, D are symmetric and commute pairwise, and K is
  their Williamson array [[X₁, X₂, X₃, X₄], [−X₂, X₁, −X₄, X₃], [−X₃, X₄, X₁, −X₂], [−X₄, −X₃, X₂, X₁]], relabelled so
  that the half swap of each quadrant becomes j ↦ 1 − j on ℤ_254.
- T-matrices of order 3 (first rows, from base sequences BS(2, 1)): `+00`, `0+0`, `00+`, `000`.

## Search

None: the matrix is determined by q, the near-Williamson matrices and the T-matrices. `python3 search/miyamoto_blocks.py
certify 509 3 nw:near-williamson/NW127.txt` writes the certificate.

## Status

**First explicit construction** ([MIYAMOTO_BLOCKS.md](../MIYAMOTO_BLOCKS.md#near-williamson-matrices)). Cati–Pasechnik's
Table 4 lists n = 1527 as open. Đoković (arXiv:1008.2043, Proposition 3.3; Int. Math. Forum 6 (2011) 3061–3067)
claims it with Williamson-type matrices of order 509, which he obtains from his good matrices of order 127 through
Seberry–Yamada's Corollary 8.8(i). The printed proof of that corollary does not give Williamson-type matrices
([NOTES.md](../../NOTES.md) §7), and [FRONTIER.md](../../FRONTIER.md) listed the order as unresolved. Good matrices of
order 127 also satisfy the hypotheses used here (circulant, B, C, D symmetric), so Đoković's matrices would serve as
well as the near-Williamson matrices of the 2026 database. Miyamoto's 1991 paper (J. Combin. Theory Ser. A 57) may record
the order too: by the transcription in the JD-Jones-ASES Hadamard-M note, its list (2) of Williamson-type orders
contains 127, which would give 509 and so 1527 = 3 · 509 in its list (3); the same note finds that entry unsupported by
Miyamoto's Corollary 4. We could not read the paper itself (zbMATH review Zbl 0733.05017), so this is unverified.

- Paley's constructions do not apply: 4n − 1 = 6107 = 31 · 197 and 2n − 1 = 3053 = 43 · 71 are not prime powers.
- Miyamoto's theorem does not apply: it needs a symmetric conference matrix of order n + 1, which cannot exist for n ≡ 3
  (mod 4).
- Yamada–Kiyasu does not apply: n − 2 = 1525 = 5² · 61 is not a prime power ≡ 1 (mod 4).
- Cooper–Wallis needs Williamson-type matrices of order 509 or T-matrices of order 509; neither is known to us.
- The Paley-block route does not reach it: 4s − 1 = 2035 and 6107 are not prime powers ≡ 3 (mod 8) for the other
  divisors s of 1527, and the remaining ones need T-matrices of order 1527 (s = 1) or 509 (s = 3), which are not known
  to us.
- The Miyamoto-block route with a negaperiodic or periodic Golay pair of length 254 does not apply: no negaperiodic pair
  of that length is known to us (Balonin–Đoković's lists (25) and (31) cover odd orders below 120 only, and none of their
  families gives length 254: 253 = 11 · 23 and 507 = 3 · 13² are not prime powers), and a periodic one cannot exist (508 = 2² · 127 is not a sum of two squares).
- SageMath 10.8 (passagemath 10.8.13, October 2026) does not construct it: `hadamard_matrix(6108, existence=True)`
  returns Unknown.

We found no valid construction of order 6108 before this one. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the near-Williamson matrices, the T-matrices and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py 6108
```

SHA-256 of `H6108.csv`: `98a99283752d0b9c431f84ea315675b86a8fd1592a863c14ea26f35c969e2d71`

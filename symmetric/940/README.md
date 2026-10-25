# Symmetric Hadamard matrix of order 940

940 = 4 · 235, with 235 = 5 · 47; 235 = q + 2 with q = 233 a prime ≡ 1 (mod 4) and (q + 1)/2 = 117.

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is Seberry and Balonin's Corollary 1 (Australas. J. Combin. 69
(2017), Theorem 3 and Corollary 1, with the order corrected by Balonin and Đoković, arXiv:1803.08787), as for
[symmetric/1804](../1804): a symmetric Hadamard matrix of order 4(q + 2) from the core of a symmetric conference matrix
of order (q + 3)/2, here 118.

- **The conference matrix of order 118.** A symmetric conference matrix C of order 118 (zero diagonal, entries ±1
  elsewhere, C² = 117 I), found here (round 44, [search/NOTES.md](../../search/NOTES.md) §95) on the point set
  ℤ₁₃ × ℤ₃² plus one point: its 9 × 9 blocks are hyperplane matrices of AG(2, 3) or their complements, and the block
  pattern is an exact cover of the 78 pairs of points of ℤ₁₃ by 13 "K₄ components" (the lines of PG(2, 3)), with
  Paley signs on ℤ₁₃. The matrix is `search/data/conference118_5fe566972f21790b.txt` (normalised: first row and column
  +, zero diagonal); its core Q of order 117 is in the certificate. Its signed automorphism group is {±I}.
- **Yamada's circulants.** M and N of order 117 with M² + N² = 233 I (Yamada's Lemma 3 over GF(233²), r = 3 a
  non-square, ξ as in the certificate).
- **Williamson-type matrices.** Seberry–Yamada's Lemma 26 with U = (I, 0, Q, Q) and V = (M, N, I, I) gives symmetric,
  pairwise amicable W₁, …, W₄ of order 235 with Σ W_k² = 940 I and W₃ = W₄ (their Corollary 27).
- **The matrix.** H = [[W₁, W₃, W₃, W₂], [W₃, W₂, −W₁, −W₃], [W₃, −W₁, −W₂, W₃], [W₂, −W₃, W₃, −W₁]], Seberry and
  Balonin's propus array.

## Search

The conference matrix came from an exact cover (search/NOTES.md §95): the point pairs of ℤ₁₃ are covered by components
of a cubic graph on 104 flag vertices, each component a K₄ with signs, ±1 diagonal and cube-root phases chosen from 48
solvable classes; the first cover appears after 13,589 nodes (0.09 s). The rest is deterministic.

| Code | Command | Randomness and caps | Result |
|---|---|---|---|
| `search/sym_propus_sb.py`, `verify/verify.py` construction `symmetric-propus-sb` | `python3 -c "import numpy as np; np.save('CONF.npy', np.array([[{'0': 0, '+': 1, '-': -1}[c] for c in r] for r in open('search/data/conference118_5fe566972f21790b.txt').read().split()]))"`, then `python3 search/sym_propus_sb.py certify 235 CONF.npy` | none | `certificate.json`, matrix SHA-256 `dad8b704…` |

## Status

**Concurrent.** The Đoković–Seberry census (2025-05-28) lists 235 among its undecided values, and we found no
construction of order 940 by anyone else.

- The census's Seberry–Balonin series needs (q + 1)/2 = 117 to be a prime power or the order of the core of a symmetric
  conference matrix of order 118. Brouwer's table of strongly regular graphs (read on 17 October 2026) marks the
  conference graph SRG(117, 58, 28, 29) as existing, citing J. Y. Z. Ho (email, October 2026); that matrix is not
  public, and we do not know whether Ho derived this order. Ours was found independently. Since Ho's graph gives this
  order through the published Corollary 1, we label it Concurrent and claim no priority, as for
  [hadamard/2284](../../hadamard/2284/README.md#status). [hadamard/MIYAMOTO_BLOCKS.md](../../hadamard/MIYAMOTO_BLOCKS.md)
  wrote this chain down before the conference matrix was known to us.
- No other route that we know: the propus searches of [search/NOTES.md](../../search/NOTES.md) §17 and §29 and the
  plug-in searches of §87 and §94 found nothing for v = 235. 235 = 5 · 47 is not a prime power, so the constructions of
  order 4q of this directory ([PROPUS_MIYAMOTO.md](../PROPUS_MIYAMOTO.md), [PROPUS_WEIGHING.md](../PROPUS_WEIGHING.md),
  [PERIODIC_GOLAY.md](../PERIODIC_GOLAY.md)) do not apply; Yamada–Kiyasu ([YAMADA_KIYASU.md](../YAMADA_KIYASU.md))
  needs q = 233 ≡ 5 (mod 8), but 233 ≡ 1; and 940 is not of the form 2n(n − 1), m(m − 1) or m(m + 3).
- Ordinary Hadamard matrices of order 940 are classical.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, r, ξ, the core Q of order 117 (its strict upper triangle, base64, 1 for +1) with its SHA-256, and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/940
```

`verify/verify.py` (construction `symmetric-propus-sb`) unpacks Q and checks its SHA-256, that it is symmetric with zero
diagonal, Qe = 0 and Q² = 117I − J; builds Yamada's circulants and checks M² + N² = 233I; builds the W_k and checks that
they are ±1, symmetric, pairwise amicable, with Σ W_k² = 940 I and W₃ = W₄; and checks H = Hᵀ, H Hᵀ = 940 I and the
SHA-256 of the canonical CSV.

SHA-256 of `H940.csv`: `dad8b7046baf706b070c2706de9c65c2725604119dd026d008f47e88f0c750d5`

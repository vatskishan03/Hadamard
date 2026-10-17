# Symmetric Hadamard matrix of order 1804

1804 = 4 · 451, with 451 = 11 · 41; 451 = q + 2 with q = 449 a prime ≡ 1 (mod 4) and (q + 1)/2 = 225.

## Construction

A Hadamard matrix H is symmetric if H = Hᵀ. This one is Seberry and Balonin's Corollary 1 (Australas. J. Combin. 69
(2017), Theorem 3 and Corollary 1; arXiv:1512.01732), the Seberry–Balonin series of the Đoković–Seberry census: for a
prime power q ≡ 1 (mod 4) such that (q + 1)/2 is a prime power or the order of the core of a symmetric conference
matrix, there is a symmetric Hadamard matrix of order 4(q + 2). Here (q + 1)/2 is not a prime power, and the core comes
from

- **The core.** a symmetric conference matrix of order 226 = 15² + 1 from a skew-Hadamard matrix of order 16 (Brouwer's
  table of strongly regular graphs lists the conference graph SRG(225, 112, 55, 56) as "skewhad²"): with the skew core S
  of order 15 of Paley's skew matrix of order 8 doubled (Sᵀ = −S, SJ = 0, SSᵀ = 15I − J), Q = S ⊗ S + I ⊗ J − J ⊗ I is
  symmetric with zero diagonal, Qe = 0 and Q² = 225I − J.
- **Yamada's circulants.** M and N of order 225, with first rows β and α of Yamada's Lemma 3 over GF(449²) =
  GF(449)[w]/(w² − 3) and ξ = 6 + w primitive (as in [hadamard/15292](../../hadamard/15292)): symmetric, M of zero
  diagonal, M² + N² = 449I.
- **Williamson-type matrices.** Seberry–Yamada's Lemma 26 with U = (I, 0, Q, Q) and V = (M, N, I, I) gives W₁ = [[1,
  −e], [−eᵀ, S₁]] and W_k = [[1, e], [eᵀ, S_k]] (k > 1), S_k = U_k ⊗ J₂ + V_k ⊗ [[1, −1], [−1, 1]], of order 451:
  symmetric, pairwise amicable, Σ W_k² = 1804 I, and W₃ = W₄. This is also Seberry–Yamada's Corollary 31 (J. Combin.
  Math. Combin. Comput. 7 (1990)) with a conference core in place of Paley's.
- **The matrix.** H = [[W₁, W₃, W₃, W₂], [W₃, W₂, −W₁, −W₃], [W₃, −W₁, −W₂, W₃], [W₂, −W₃, W₃, −W₁]], Seberry and
  Balonin's propus array.

## Search

None. `python3 search/sym_propus_sb.py skewsq 16 conf226.npy` writes the conference matrix, and `python3
search/sym_propus_sb.py certify 451 CONF.npy` the certificate.

## Status

**Known**: the order follows from Seberry and Balonin's Corollary 1 (2015/2017) with a conference matrix published long
before, but the Đoković–Seberry census (2025-05-28) omits it from its Seberry–Balonin series and lists 451 among its 595
undecided values. We found no explicit matrix.

- The census's Seberry–Balonin series is "v = q + 2, q ≡ 1 (mod 4), and (q + 1)/2 a prime power or the order of the core
  of a symmetric conference matrix"; q = 449 and (q + 1)/2 = 225 satisfy it. Seberry and Balonin name q = 89 (core of
  order 45, Mathon's matrix of order 46) as an instance of the second kind. Among the symmetric conference matrices we
  know whose core order is not a prime power (Mathon's orders 46, 442 and 1574, and 226 from a skew-Hadamard matrix of
  order 16), only those of orders 46, 226 and 442 give v < 2500: v = 91, 451 and 883.
- Ordinary Hadamard matrices of order 1804 are classical, and a skew-Hadamard matrix of this order is in
  [skew/1804](../../skew/1804).
- Seberry and Balonin print the order of Corollary 1 as 4(2q + 1), with q in place of (q + 1)/2; their examples (11 and
  83) and the census read it as 4(q + 2).

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, r, ξ, the core Q of order 225 (its strict upper triangle, base64, 1 for +1) with its SHA-256, and the SHA-256 of the matrix |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py symmetric/1804
```

`verify/verify.py` (construction `symmetric-propus-sb`) unpacks Q and checks its SHA-256, that it is symmetric with zero
diagonal, Qe = 0 and Q² = 225I − J; builds Yamada's circulants (checking that r is a non-square and ξ primitive) and
checks M² + N² = 449I; builds the W_k and checks that they are ±1, symmetric, pairwise amicable, with Σ W_k² = 1804 I
and W₃ = W₄; and checks H = Hᵀ, H Hᵀ = 1804 I and the SHA-256 of the canonical CSV.

SHA-256 of `H1804.csv`: `ad1f51fd81c483148f332a78267737753d1cbd0242c8bc68b47d1b84a4402020`

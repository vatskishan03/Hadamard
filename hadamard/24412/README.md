# Order 24412

24412 = 4 · 6103, with 6103 = 17 · 359.

## Construction

**Yamada's theorem.** M. Yamada, *Some new series of Hadamard matrices*, J. Austral. Math. Soc. (Ser. A)
46 (1989) 371–383, Theorem 4, which Yamada credits to Z. Kiyasu: if q ≡ 5 (mod 8) is a prime power and a
skew-Hadamard matrix of order (q + 3)/2 exists, then a Hadamard matrix of order 4(q + 2) exists.

Here q = 6101, which is prime and ≡ 5 (mod 8), so (q + 3)/2 = 3052 and 4(q + 2) = 24412.

**The skew input** is this repository's skew-Hadamard matrix of order 3052 = 4 · 763
([skew/3052](../../skew/3052)), which is new. Before it, the smallest known skew order with odd part 763 was
2⁴ · 763 = 12208 (Cati–Pasechnik Table 5: 763(4)), so a skew matrix of order 3052 was not available.

**The array.** The construction follows Yamada's proof:

- the circulants α, β of order (q + 1)/2 = 3051 from the quadratic character of GF(q) on traces of powers
  of a primitive element of GF(q²) (Yamada's relative Gauss sums, Lemma 3: α αᵀ + β βᵀ = qI);
- the normalised skew matrix [[1, eᵀ], [−e, S + I]];
- the blocks A = P ⊗ β + J₂ ⊗ I, B = P ⊗ α, C = D = J₂ ⊗ S + P ⊗ I;
- the bordered quaternion-type array of Yamada's Theorem 1.

`build_yamada` in `verify/verify.py` implements this, as for orders [24124](../24124) and
[51580](../51580). The certificate gives the inputs:
- q = 6101;
- the non-square r = 2, which defines GF(q²) = GF(q)[w]/(w² − r);
- the primitive element ξ = 2 + w;
- the nested certificate of the skew matrix of order 3052.

## Status

**Probably new.** The checks:

- **Handbook of Combinatorial Designs (2007), Table 1.53, lists n = 6103 with t = 3.** So the smallest
  known Hadamard matrix with odd part 6103 had order 8 · 6103. Paley I gives that one, since
  8 · 6103 − 1 = 48823 is prime. We read the table from an unofficial copy and checked the parse against
  known entries (6031 → 3, 5437 → 8, and every value Paley reaches shows no entry). Treat the claim with
  that caveat.
- **Later updates:** Đoković's updates (arXiv:0912.5091, arXiv:1008.2043) do not include 6103.
  Cati–Pasechnik Table 4 stops at n = 2999.
- **Other routes** do not apply:
  - Paley I and II: 24411 = 3 · 79 · 103 and 12205 = 5 · 2441.
  - Miyamoto: 6103 is not a prime power and is ≡ 3 (mod 4).
  - Yamada's other two theorems: these need q ≡ 1 (mod 8).
  - Williamson-type families of order (q + 1)/2: 717 and 12205 are not prime powers.
  - Cooper–Wallis with 6103 = 17 · 359: this needs Williamson-type matrices of order 359 or T-matrices of
    order 359. Neither is known.
  - Kronecker and conference-matrix products: no admissible split.
- **2026 deposits:** none mentions 24412 or 6103. We checked Alpöge et al., Schneider, the JD-Jones-ASES
  repositories Hadamard-T/M/B, Schumer, Karoui and Piza.

This is a "found nothing" result, not a proof that no earlier construction exists. The literature audit
was done by a research agent on 2026-10-12.

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, r, ξ and the nested certificate of the skew matrix of order 3052 |

The matrix is not stored, because its CSV encoding takes about 1.2 GB.

SHA-256 of `H24412.csv`: `2f9eb90bc7fd0ab3d76f8274af24ad0e4d4c48000225e126df0d5388b52be3c7`.
It was built and checked by `python3 search/make_yamada.py 3052 0` in 4 minutes on 4 shared cores.

```
python3 verify/verify.py 24412
```

rebuilds the skew input from its certificate and applies Yamada's construction. It then checks
H Hᵀ = 24412·I exactly, with float32 blocks, and compares the SHA-256 of the canonical CSV.

A separate run of this command (2026-10-12, 4 minutes) reported `24412 ok  H H^T = 24412 I` with the
digest above.

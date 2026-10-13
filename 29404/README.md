# Order 29404

29404 = 4 · 7351, with 7351 prime and 7351 ≡ 3 (mod 4).

## Construction

**Yamada's theorem.** M. Yamada, *Some new series of Hadamard matrices*, J. Austral. Math. Soc. (Ser. A)
46 (1989) 371–383, Theorem 4, which Yamada credits to Z. Kiyasu: if q ≡ 5 (mod 8) is a prime power and a
skew-Hadamard matrix of order (q + 3)/2 exists, then a Hadamard matrix of order 4(q + 2) exists.

Here q = 7349, which is prime and ≡ 5 (mod 8), so (q + 3)/2 = 3676 and 4(q + 2) = 29404. (7349 and 7351 are
twin primes.)

**The skew input** is this repository's skew-Hadamard matrix of order 3676 = 4 · 919 ([skew/3676](../skew/3676)).
That order was first published by Schumer (Zenodo, doi:10.5281/zenodo.22820997, 2026-09-17); our matrix is an
independent construction with the same spin structure. Cati–Pasechnik list 919(3) in both Table 4 and Table 5,
so no skew matrix of order 3676 was available before September 2026.

**The array** is the one of [order 24412](../24412): the circulants α, β of order (q + 1)/2 = 3675 from the
quadratic character of GF(q) on traces of powers of a primitive element of GF(q²), the normalised skew matrix
[[1, eᵀ], [−e, S + I]], the blocks A = P ⊗ β + J₂ ⊗ I, B = P ⊗ α, C = D = J₂ ⊗ S + P ⊗ I, and Yamada's bordered
quaternion-type array. `build_yamada` in `verify/verify.py` implements it. The certificate gives q = 7349, the
non-square r = 2 (GF(q²) = GF(q)[w]/(w² − 2)), the primitive element ξ = 12 + w and the nested certificate of the
skew matrix of order 3676.

## Status

First explicit construction we know of, with one caveat.

- **Caveat: Handbook of Combinatorial Designs (2007), Table 1.53, marks n = 7351 with "." (t = 2), that is, as
  known.** We found no construction behind that entry:
  - Paley I and II fail: 29403 = 3⁵ · 11² and 14701 = 61 · 241.
  - Miyamoto does not apply (7351 ≡ 3 mod 4).
  - Yamada's theorem with q = 7349 needed a skew matrix of order 3676, which did not exist before 2026.
  - No Williamson-type family in the Seberry–Yamada (1992) legend reaches 7351, and Spence's 4(q² + q + 1)
    does not fit.
  - Classical theorems give 8 · 7351 (Miyamoto's 8q theorem, since 2 · 7351 − 3 = 14699 is prime, and Spence's
    twin-prime-power theorem), not 4 · 7351.
  - Yamada's own 1989 list of new orders from this theorem (n < 10000) does not contain 7351.
- **Why the entry looks like an error.** Between 4000 and 9999, every twin-prime n ≡ 3 (mod 4) without a Paley
  route is marked "." (39 of 39), although only 6 of them have a construction we could identify; non-twin primes in
  the same class are marked "." in 3 of 195 cases. Below 4000 the "." twin entries all have real constructions. The
  table seems to have applied a twin-prime rule that gives 4n only for n ≡ 1 (mod 4).
- **Later updates and deposits.** Đoković's updates (arXiv:0912.5091, arXiv:1008.2043) and Cati–Pasechnik
  (n ≤ 2999) do not cover it; no 2024–2026 deposit we checked (GitHub, Zenodo, arXiv, the 2026 constructions of
  Alpöge et al., Schneider, JD-Jones-ASES, Schumer) mentions 29404 or 7351. SageMath (passagemath 10.8.13) returns
  Unknown. The order has followed from Schumer's skew matrix with Yamada's theorem since 2026-09-17, but nobody
  states it.

We read Table 1.53 from an unofficial copy and checked the parse against the whole table (Paley entries, Đoković's
corrections and his count of 1006 non-dot entries). This is a "found nothing" result, not a proof that no earlier
construction exists; the audit was done by a research agent on 2026-10-13.

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, r, ξ and the nested certificate of the skew matrix of order 3676 |

The matrix is not stored (its CSV encoding takes about 1.7 GB).

SHA-256 of `H29404.csv`: `f8edf5363a448cc92377e0b3293199bab28e94e8d541f7b1335ba87fc1f8e3ef`.
It was built and checked (H Hᵀ = 29404 I) by `python3 attack/make_yamada.py 3676 0`; `python3 verify/verify.py 29404`
rebuilds and checks it from the certificate in 8 minutes on 4 shared cores (2026-10-13).

```
python3 verify/verify.py 29404
```

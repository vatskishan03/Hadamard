# Order 51580

51580 = 4 · 12895, with 12895 = 5 · 2579.

## Construction

**Yamada's theorem.** M. Yamada, *Some new series of Hadamard matrices*, J. Austral. Math. Soc. (Ser. A)
46 (1989) 371–383, Theorem 4, which Yamada credits to Z. Kiyasu: if q ≡ 5 (mod 8) is a prime power and a
skew-Hadamard matrix of order (q + 3)/2 exists, then a Hadamard matrix of order 4(q + 2) exists.

Here q = 12893, which is prime and ≡ 5 (mod 8), so (q + 3)/2 = 6448 and 4(q + 2) = 51580.

**The skew input.** The skew-Hadamard matrix of order 6448 = 16 · 403 is built in two steps:

- this repository's skew-Hadamard matrix of order 1612 = 4 · 403 ([skew/1612](../skew/1612)), which is new;
- two doublings [[S, S], [−Sᵀ, Sᵀ]] (Seberry), giving orders 3224 and 6448.

Before 1612, the smallest known skew order with odd part 403 was 2⁵ · 403 = 12896 (Cati–Pasechnik
Table 5: 403(5)), so a skew matrix of order 6448 was not available.

**The array.** The construction follows Yamada's proof:

- the circulants α, β of order (q + 1)/2 = 6447 from the quadratic character of GF(q) on traces of powers
  of a primitive element of GF(q²) (Yamada's relative Gauss sums, Lemma 3: α αᵀ + β βᵀ = qI);
- the normalised skew matrix [[1, eᵀ], [−e, S + I]];
- the blocks A = P ⊗ β + J₂ ⊗ I, B = P ⊗ α, C = D = J₂ ⊗ S + P ⊗ I;
- the bordered quaternion-type array of Yamada's Theorem 1.

`build_yamada` in `verify/verify.py` (from main) implements this. Its inputs are q = 12893, the non-square
r = 2 that defines GF(q²) = GF(q)[w]/(w² − r), the primitive element ξ = 3 + w, and the nested
certificate of the skew matrix.

## Status

**Probably new.** The checks:

- 12895 > 2999 lies beyond Cati–Pasechnik Table 4, and 12895 > 10000 is beyond Table 1.53 of the
  Handbook of Combinatorial Designs (2007), which covers odd n < 10000.
- **Paley:** 51579 and 25789 are composite, so neither Paley construction applies.
- **Williamson-type routes:** in Seberry–Yamada's 1992 survey, Table A.2 gives n = 2579 only at order
  8 · 2579. No family in its legend produces Williamson-type matrices of order 2579 or 12895. Base
  sequences, Golay sums and Yang numbers do not reach T-sequence lengths 2579 or 12895.
- **6m + 1 T-matrices:** 12895 is not prime, so that method cannot reach it.
- **Kiyasu–Yamada with older inputs:** this route needs a skew matrix of order 6448 = 16 · 403, which
  needs the skew order 1612 of this repository.
- **2026 deposits:** none mentions the order. We checked Alpöge et al., Schneider, JD-Jones-ASES
  Hadamard-T/M/B, Schumer, Karoui, Piza, and a scan of about 175 recent Zenodo records.

This is a "found nothing" result, not a proof that no earlier construction exists, and it carries that
caveat. The literature audit was done by research agents on 2026-10-12 and is summarised in
attack/NOTES.md §13.

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, r, ξ and the nested certificates: skew 1612, then two doublings |

The matrix is not stored. Its CSV encoding takes about 7 GB.

SHA-256 of `H51580.csv`: `de5fd6da01bf62b744026b8e92169e75eb8f36176985e578fa5150439de71cbd`. Built and checked by `python3 attack/make_yamada.py 1612 2` in 44 minutes on 4 shared cores, about 3 GB of memory.

```
python3 verify/verify.py 51580
```

rebuilds the skew input from `skew/1612/certificate.json`, doubles it twice, applies Yamada's
construction, checks H Hᵀ = 51580·I exactly (float32 blocks, as for main's Miyamoto orders) and compares
the SHA-256 of the canonical CSV. A separate run of this command (2026-10-12, 40 minutes, about 3.5 GB of memory) reported
`51580 ok  H H^T = 51580 I` with the digest above.

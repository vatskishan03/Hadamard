# Skew-Hadamard matrix of order 13072

13072 = 16 · 817 = 172 · 76, with 817 = 19 · 43.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I. This one is J. Wallis's product (J. Combin. Theory
Ser. A 21 (1976); Seberry–Yamada (1992), SXI) of a skew-Hadamard matrix of order 172 = 4 · 43 with the amicable pair of
order 76 of [skew/5700](../5700).

- **A skew-Hadamard matrix of order 172.** Four circulants of order 43 in the Goethals–Seidel array: A₀ has −1 on the
  nonzero squares mod 43 (Paley; skew type), and A₁, A₂, A₃ are unions of orbits of
  H = ⟨4⟩ = {1, 4, 11, 16, 21, 35, 41}, with row sums 1, −1, 13 and PAF_A₁ + PAF_A₂ + PAF_A₃ = 1 off zero, so the
  four autocorrelations sum to 172δ. This is the class of difference families of Đoković and Kotsireas
  (arXiv:1511.08734), who list v = 43 among the solved parameter sets. Skew-Hadamard matrices of order 172 were first
  given by Đoković (J. Combin. Theory Ser. A 61 (1992)), and any of them would do.
- **The amicable pair of order 76.** A = I + PN (skew-Hadamard) and M = P + N (symmetric Hadamard), with A Mᵀ = M Aᵀ,
  from the skew-propus quadruple of length 19 of [skew/5700](../5700) (step 1 of its README).
- **The product.** With I + S the matrix of order 172 (S skew, S Sᵀ = 171 I),

      K = I ⊗ A + S ⊗ M        (order 172 · 76).

  K + Kᵀ = I ⊗ (A + Aᵀ) + (S + Sᵀ) ⊗ M = 2I, and since A Mᵀ = M Aᵀ,

      K Kᵀ = I ⊗ A Aᵀ + (S + Sᵀ) ⊗ A Mᵀ + S Sᵀ ⊗ M Mᵀ = 76 I + 171 · 76 I = 13072 I.

## Search

`python3 search/skew_wallis.py paley3 43 7` enumerates every union of the seven orbits of ⟨4⟩ on ℤ₄₃ (orbit {0}
included), tabulates them by autocorrelation and joins three; it finishes in about a second, and the certificate uses
its first solution. The amicable pair is that of [skew/5700](../5700).

**Reproduction.**

| | |
|---|---|
| Code | `search/skew_wallis.py` (`paley3`, `certify`) and `verify/verify.py` (`build_skew_amicable_product`) |
| Command | `python3 search/skew_wallis.py certify 13072` writes the certificate; `python3 verify/verify.py skew/13072` rebuilds and checks it |
| Randomness and caps | none |

## Status

**New.** Cati–Pasechnik Table 5 gives 817(5): the smallest skew-Hadamard order known with odd part 817 was
32 · 817 = 26144. The Đoković–Seberry census lists 817 as undecided, so 4 · 817 = 3268 is open (our searches for it
failed, [search/NOTES.md](../../search/NOTES.md)), and so is 8 · 817 = 6536.

- The route needs the amicable pair of order 76, which we constructed in 2026 ([skew/5700](../5700#status)). Every
  split 13072 = h m into a skew order h and a classical amicable order m fails. m = 2 and m = 4 need the open skew
  orders 6536 and 3268. m = 8, 16, 152, 304, 344 and 688 leave h ≢ 0 (mod 4). The smallest classical amicable orders
  with odd parts 19 and 43 are 152 and 344 (Seberry–Yamada (1992), Table 6.1; Paley with q = 151 and 343).
- Our closure of the classical skew theorems also stops at 2⁵ · 817, as Table 5 does. It covers Paley, Szekeres,
  Whiteman, Wallis–Whiteman, Spence, Turyn, Seberry's theorems with classical amicable orders and orthogonal designs,
  Wallis's product and doubling ([NOTES.md](../../NOTES.md) §28). It does not cover Seberry–Yamada's SX, or Seberry's
  theorem with m ≥ 2, for which we know no classical family of orthogonal designs.
- In SageMath's development version of 2026-10-17, `skew_hadamard_matrix(13072, existence=True)` is `Unknown`; it is
  `True` for 26144. Its product branch uses only amicable orders d with d − 1 a prime power.
- Of the ingredients, Wallis's theorem and the skew order 172 were published; the amicable pair of order 76 is new.

Doubling (Seberry) gives skew-Hadamard matrices of all orders 2ᵏ · 13072. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | the factor of order 172 ("factor": v, multiplier, orbit representatives of the four blocks, row sums), the skew-propus quadruple of the pair ("amicable"), and the SHA-256 of the matrix |
| `sequences.txt` | A₀, A₁, A₂, A₃ (length 43), then the skew-propus quadruple A3, A1, A0 of [skew/5700](../5700) (length 19), one per line |

The matrix itself is not stored (about 420 MB as CSV); `verify/verify.py` rebuilds it in about a minute.

## Verification

```
python3 verify/verify.py skew/13072
python3 verify/verify.py skew/13072 --write && ./check_matrix < skew/13072/H13072.csv
```

`verify/verify.py` (construction `skew-amicable-product`) rebuilds the Goethals–Seidel matrix of order 172 and checks
that it is skew-Hadamard. It rebuilds the amicable pair from the skew-propus quadruple and checks A + Aᵀ = 2I,
M = Mᵀ, A Aᵀ = M M = 76 I and A Mᵀ = M Aᵀ. Then it assembles K and checks K Kᵀ = 13072 I, K + Kᵀ = 2I and the SHA-256
of the canonical CSV.

SHA-256 of `H13072.csv`: `c825ece90f9b77b846546eba902305b7f69051b9300325d526a40b5e8075f70a`

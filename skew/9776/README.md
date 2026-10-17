# Skew-Hadamard matrix of order 9776

9776 = 16 · 611 = 188 · 52, with 611 = 13 · 47.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I. This one is J. Wallis's product (J. Combin. Theory
Ser. A 21 (1976); Seberry–Yamada (1992), SXI) of Đoković's skew-Hadamard matrix of order 188 = 4 · 47 with the
amicable pair of order 52 of [skew/2652](../2652).

- **A skew-Hadamard matrix of order 188.** Đoković, "Skew-Hadamard matrices of orders 188 and 388 exist",
  arXiv:0704.0640 (Int. Math. Forum 3 (2008)), Proposition 2.1: four circulants of order 47 in the Goethals–Seidel
  array. A₀ has −1 on the nonzero squares mod 47 (Paley; skew type), and A₁, A₂, A₃ have −1 on his sets X₁, X₂, X₃
  (a 3-(47; 30, 22, 22; 39) difference family), with row sums −13, 3, 3.
- **The amicable pair of order 52.** A = I + PN (skew-Hadamard) and M = P + N (symmetric Hadamard), with A Mᵀ = M Aᵀ,
  from the skew-propus quadruple of length 13 of [skew/2652](../2652) (step 1 of its README).
- **The product.** With I + S the matrix of order 188 (S skew, S Sᵀ = 187 I),

      K = I ⊗ A + S ⊗ M        (order 188 · 52).

  K + Kᵀ = I ⊗ (A + Aᵀ) + (S + Sᵀ) ⊗ M = 2I, and since A Mᵀ = M Aᵀ,

      K Kᵀ = I ⊗ A Aᵀ + (S + Sᵀ) ⊗ A Mᵀ + S Sᵀ ⊗ M Mᵀ = 52 I + 187 · 52 I = 9776 I.

## Search

None: both factors were known before this order was built. The sets X₁, X₂, X₃ are copied from Đoković's paper
(`DJOKOVIC_47` in `search/skew_wallis.py`); `verify/verify.py` checks them. The amicable pair is that of
[skew/2652](../2652).

**Reproduction.**

| | |
|---|---|
| Code | `search/skew_wallis.py` (`certify`) and `verify/verify.py` (`build_skew_amicable_product`) |
| Command | `python3 search/skew_wallis.py certify 9776` writes the certificate; `python3 verify/verify.py skew/9776` rebuilds and checks it |
| Randomness and caps | none |

## Status

**New.** Cati–Pasechnik Table 5 gives 611(5): the smallest skew-Hadamard order known with odd part 611 was
32 · 611 = 19552. The Đoković–Seberry census lists 611 as undecided, so 4 · 611 = 2444 is open, and so is
8 · 611 = 4888.

- The route needs the amicable pair of order 52, which we constructed in 2026 ([skew/2652](../2652#status)). Every
  split 9776 = h m into a skew order h and a classical amicable order m fails. m = 2 and m = 4 need the open skew
  orders 4888 and 2444. m = 8, 16, 104, 208 and 752 leave h ≢ 0 (mod 4). The smallest classical amicable orders with
  odd parts 13 and 47 are 104 and 752 (Seberry–Yamada (1992), Table 6.1; Paley with q = 103 and 751).
- Our closure of the classical skew theorems also stops at 2⁵ · 611, as Table 5 does. It covers Paley, Szekeres,
  Whiteman, Wallis–Whiteman, Spence, Turyn, Seberry's theorems with classical amicable orders and orthogonal designs,
  Wallis's product and doubling ([NOTES.md](../../NOTES.md) §28). It does not cover Seberry–Yamada's SX, or Seberry's
  theorem with m ≥ 2, for which we know no classical family of orthogonal designs.
- In SageMath's development version of 2026-10-17, `skew_hadamard_matrix(9776, existence=True)` is `Unknown`; it is
  `True` for 19552. Its product branch uses only amicable orders d with d − 1 a prime power.
- Of the ingredients, Wallis's theorem and Đoković's skew order 188 were published; the amicable pair of order 52 is
  new.

Doubling (Seberry) gives skew-Hadamard matrices of all orders 2ᵏ · 9776. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | the factor of order 188 ("factor": v, the four sets, row sums), the skew-propus quadruple of the pair ("amicable"), and the SHA-256 of the matrix |
| `sequences.txt` | A₀, A₁, A₂, A₃ (length 47), then the skew-propus quadruple A3, A1, A0 of [skew/2652](../2652) (length 13), one per line |

The matrix itself is not stored (about 235 MB as CSV); `verify/verify.py` rebuilds it in under a minute.

## Verification

```
python3 verify/verify.py skew/9776
python3 verify/verify.py skew/9776 --write && ./check_matrix < skew/9776/H9776.csv
```

`verify/verify.py` (construction `skew-amicable-product`) rebuilds the Goethals–Seidel matrix of order 188 and checks
that it is skew-Hadamard. It rebuilds the amicable pair from the skew-propus quadruple and checks A + Aᵀ = 2I,
M = Mᵀ, A Aᵀ = M M = 52 I and A Mᵀ = M Aᵀ. Then it assembles K and checks K Kᵀ = 9776 I, K + Kᵀ = 2I and the SHA-256 of
the canonical CSV.

SHA-256 of `H9776.csv`: `9f5ab196b0219583d3b05f5c336c7cfd83f5e7ca94403e3de793ab617465498d`

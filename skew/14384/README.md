# Skew-Hadamard matrix of order 14384

14384 = 16 · 899 = 116 · 124, with 899 = 29 · 31.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I. This one is J. Wallis's product (J. Combin. Theory
Ser. A 21 (1976); Seberry–Yamada (1992), SXI) of a skew-Hadamard matrix of order 116 = 4 · 29 with an amicable pair of
order 124 = 4 · 31.

- **A skew-Hadamard matrix of order 116.** Good matrices of order 29 (Szekeres 1988, as SageMath stores them): four
  circulants A, B, C, D with A of skew type and B, C, D symmetric, row sums 1, 5, −3, 9 and
  PAF_A + PAF_B + PAF_C + PAF_D = 116δ, in the Goethals–Seidel array.
- **An amicable pair of order 124.** The skew-propus quadruple of length 31 of
  [search/NOTES.md](../../search/NOTES.md) §31 (`sym_mm1.py mpropus 31 5`: A3 symmetric, A1, A0 of skew type, constant
  on the orbits of ⟨5⟩, row sums 5, 7, 1, PAF_A3 + 2 PAF_A1 + PAF_A0 = 124δ) gives, in the Balonin array, a symmetric
  orthogonal design xP + yN of type (1, 123), and A = I + PN (skew-Hadamard) and M = P + N (symmetric Hadamard) are
  amicable, as in step 1 of [skew/2652](../2652).
- **The product.** With I + S the matrix of order 116 (S skew, S Sᵀ = 115 I),

      K = I ⊗ A + S ⊗ M        (order 116 · 124).

  K + Kᵀ = I ⊗ (A + Aᵀ) + (S + Sᵀ) ⊗ M = 2I, and since A Mᵀ = M Aᵀ,

      K Kᵀ = I ⊗ A Aᵀ + (S + Sᵀ) ⊗ A Mᵀ + S Sᵀ ⊗ M Mᵀ = 124 I + 115 · 124 I = 14384 I.

## Search

None for this order: the good matrices are copied from SageMath (`GOOD_29` in `search/skew_wallis.py`), and the
skew-propus quadruple was found by the multiplier search of [search/NOTES.md](../../search/NOTES.md) §31 (exhaustive
over the orbits of ⟨5⟩, in a fraction of a second).

**Reproduction.**

| | |
|---|---|
| Code | `search/skew_wallis.py` (`certify`), `search/sym_mm1.py` (`mpropus 31 5`) and `verify/verify.py` (`build_skew_amicable_product`) |
| Command | `python3 search/skew_wallis.py certify 14384` writes the certificate; `python3 verify/verify.py skew/14384` rebuilds and checks it |
| Randomness and caps | none |

## Status

**New.** Cati–Pasechnik Table 5 gives 899(6): the smallest skew-Hadamard order known to them with odd part 899 was
64 · 899. The Đoković–Seberry census lists 899 as undecided, so 4 · 899 = 3596 is open, and so is 8 · 899 = 7192.

- The best classical route we know gives 32 · 899: Seberry–Yamada's amicable pair of order 248 = 2²(61 + 1) (AV,
  61 = 5² + 4 · 3²) times a skew-Hadamard matrix of order 116 (a Table 5 omission found by the closure audit,
  [NOTES.md](../../NOTES.md) §28).
- The route here needs the amicable pair of order 124, which is ours (2026). Classically the smallest amicable order
  with odd part 31 is 248, and with odd part 29 it is 464 (Paley with q = 463 gives 4 · 29 · 4). Every split
  14384 = h m into a skew order h and a classical amicable order m fails. m = 2 and m = 4 need the open skew orders
  7192 and 3596. m = 8, 16, 248 and 496 leave h ≢ 0 (mod 4), and m = 464 leaves h = 31.
- In SageMath's development version of 2026-10-17, `skew_hadamard_matrix(n, existence=True)` is `Unknown` for
  n = 14384 and for 28768 = 32 · 899.
- Of the ingredients, Wallis's theorem and the good matrices of order 29 are classical; the amicable pair of order 124
  is new as far as we know (the classical amicable orders q + 1, 2(q + 1) and 2ʳ(q + 1) all fail at 124, and
  Seberry–Yamada's Table 6.1 gives t = 3 for q = 31).

Doubling (Seberry) gives skew-Hadamard matrices of all orders 2ᵏ · 14384. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | the factor of order 116 ("factor": v, the four sets, row sums), the skew-propus quadruple of length 31 ("amicable"), and the SHA-256 of the matrix |
| `sequences.txt` | A, B, C, D (length 29), then the skew-propus quadruple A3, A1, A0 (length 31), one per line |

The matrix itself is not stored (about 520 MB as CSV); `verify/verify.py` rebuilds it in about a minute.

## Verification

```
python3 verify/verify.py skew/14384
python3 verify/verify.py skew/14384 --write && ./check_matrix < skew/14384/H14384.csv
```

`verify/verify.py` (construction `skew-amicable-product`) rebuilds the Goethals–Seidel matrix of order 116 and checks
that it is skew-Hadamard. It rebuilds the amicable pair from the skew-propus quadruple and checks A + Aᵀ = 2I, M = Mᵀ, A
Aᵀ = M M = 124 I and A Mᵀ = M Aᵀ. Then it assembles K and checks K Kᵀ = 14384 I, K + Kᵀ = 2I and the SHA-256 of the
canonical CSV.

SHA-256 of `H14384.csv`: `38e67ac80946b9c3798b22275ac8ae532a7afc505eca9558db8c946d36a9a2b0`

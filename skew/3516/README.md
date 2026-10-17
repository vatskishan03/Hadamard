# Skew-Hadamard matrix of order 3516

3516 = 4 · 879, with 879 = 3 · 293. 3516 = 3(q + 1) with q = 1171, a prime ≡ 3 (mod 8).

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I. This one is the case q = 1171 of the construction in
[PALEY_SKEW.md](../PALEY_SKEW.md), which gives a skew-Hadamard matrix of order 3(q + 1) for every prime power q ≡ 3 (mod
4).

- **Paley's matrix in sixteen blocks.** Paley's skew-Hadamard matrix of order q + 1 = 1172 on PG(1, 1171), signed and
  listed along the four orbits, of length s = 293, of g = (x ↦ −1/(x + 4)) as in
  [PALEY_BLOCKS.md](../../hadamard/PALEY_BLOCKS.md) (sign c = −1), is a 4 × 4 array of circulant s × s blocks N_kj.
- **The pairing.** After the base point of orbit 3 is moved by 231 (`shifts`), the blocks of column 0 in rows 2 and 3
  satisfy N_30 = N_20ᵀ. The pairing comes from an anti-automorphism of Paley's matrix that fixes two of the orbits and
  swaps the other two.
- **T-matrices of order 3.** T_0 = I, T_2 = Z, T_3 = −Z² and T_1 = 0, with Z the cyclic shift (`T`, one row per block
  row). The matrices A_j = Σ_k T_k ⊗ N_kj are ±1 matrices on ℤ₈₇₉ = ℤ₃ × ℤ₂₉₃, Σ_j A_j A_jᵀ = 3516 I, and A_0 = I ⊗ N_00
  + Z ⊗ N_20 − Z² ⊗ N_20ᵀ is skew-type.
- **The matrix.** The Goethals–Seidel array with A_0 first and R the back-diagonal identity is a skew-Hadamard matrix of
  order 3516. `sequences.txt` lists the four A_j as sequences on ℤ₈₇₉, A_0 first.

## Search

None: the matrix is determined by q. `python3 search/paley_blocks.py skew 1171` writes the certificate, and `python3
search/skew_3q.py 1171` builds a matrix of the same family from the closed formula of [PALEY_SKEW.md](../PALEY_SKEW.md)
and checks it.

## Status

**First explicit construction** ([PALEY_SKEW.md](../PALEY_SKEW.md#status)). The Đoković–Seberry census lists 879 as
undecided, but two older tables imply that a skew-Hadamard matrix of order 4 · 879 = 3516 is known, and we found no
construction behind either entry, as for [2652](../2652).

- The census lists 879 among its 737 undecided v, in none of its Paley, Szekeres, Szekeres–Whiteman or Spence series and
  not among its 41 sporadic values.
- Seberry–Yamada (1992), Table 7.1, leaves 879 blank, and "a blank represents 2"; Cati–Pasechnik Table 5 (which omits m
  = 2, and whose footnote lists its own changes to Seberry–Yamada's 2020 table) has no entry for 879. Both tables treat
  663 and 783 the same way. No construction is cited, and none of the classical families gives 3516: 4 · 879 − 1 = 3515
  = 5 · 19 · 37 and 2 · 879 − 1 = 1757 = 7 · 251 are not prime powers. For all three values 2v − 1 ≡ 5 (mod 8) is
  composite, which suggests (but does not prove) a misapplied Szekeres entry 2(q + 1).

Doubling (Seberry) gives skew-Hadamard matrices of all orders 2ᵏ · 3516. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power and c (the blocks), shifts and columns (the pairing), the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences, one per line |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py skew/3516
```

SHA-256 of `H3516.csv`: `f0d83942814f7150778b67e0c1b6e0d28ecc6bb7d9ec5a04a075533f231584d8`

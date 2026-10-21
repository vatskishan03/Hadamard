# Skew-Hadamard matrix of order 2724

2724 = 4 · 681, with 681 = 3 · 227. 2724 = 3(q + 1) with q = 907, a prime ≡ 3 (mod 8).

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I. This one is the case q = 907 of the construction in
[PALEY_SKEW.md](../PALEY_SKEW.md), which gives a skew-Hadamard matrix of order 3(q + 1) for every prime power q ≡ 3 (mod
4).

- **Paley's matrix in sixteen blocks.** Paley's skew-Hadamard matrix of order q + 1 = 908 on PG(1, 907), signed and
  listed along the four orbits, of length s = 227, of g = (x ↦ −1/(x + 3))² as in
  [PALEY_BLOCKS.md](../../hadamard/PALEY_BLOCKS.md) (sign c = 1), is a 4 × 4 array of circulant s × s blocks N_kj.
- **The pairing.** After the base point of orbit 3 is moved by 76 (`shifts`), the blocks of column 0 in rows 1 and 3
  satisfy N_30 = N_10ᵀ. The pairing comes from an anti-automorphism of Paley's matrix that fixes two of the orbits and
  swaps the other two.
- **T-matrices of order 3.** T_0 = I, T_1 = Z, T_3 = −Z² and T_2 = 0, with Z the cyclic shift (`T`, one row per block
  row). The matrices A_j = Σ_k T_k ⊗ N_kj are ±1 matrices on ℤ₆₈₁ = ℤ₃ × ℤ₂₂₇, Σ_j A_j A_jᵀ = 2724 I, and A_0 = I ⊗ N_00
  + Z ⊗ N_10 − Z² ⊗ N_10ᵀ is skew-type.
- **The matrix.** The Goethals–Seidel array with A_0 first and R the back-diagonal identity is a skew-Hadamard matrix of
  order 2724. `sequences.txt` lists the four A_j as sequences on ℤ₆₈₁, A_0 first.

## Search

None: the matrix is determined by q. `python3 search/paley_blocks.py skew 907` writes the certificate, and `python3
search/skew_3q.py 907` builds a matrix of the same family from the closed formula of [PALEY_SKEW.md](../PALEY_SKEW.md)
and checks it.

## Status

**New, but obtainable classically** as a skew-Hadamard matrix ([PALEY_SKEW.md](../PALEY_SKEW.md#status)): Yamada's
sixteen-block form of Paley's matrix with its transposed pair (1991) and Seberry and Yamada's multiplication by
T-matrices (1990) were published, and we found them combined nowhere.

- The Đoković–Seberry census lists 681 among its 737 undecided v, in none of its Paley, Szekeres, Szekeres–Whiteman or
  Spence series and not among its 41 sporadic values.
- Cati–Pasechnik Table 5 gives 681(4): the smallest known skew order with odd part 681 was 16 · 681 = 10896.
  Seberry–Yamada's Table 7.1 (1992) agrees.

Doubling (Seberry) gives skew-Hadamard matrices of all orders 2ᵏ · 2724. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power and c (the blocks), shifts and columns (the pairing), the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences, one per line |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py skew/2724
```

SHA-256 of `H2724.csv`: `c856ff1d10bc2695b0224604dc416246bb56d943fa1f6f667049c027992d1e57`

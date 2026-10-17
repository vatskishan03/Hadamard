# Skew-Hadamard matrix of order 6612

6612 = 4 · 1653, with 1653 = 3 · 19 · 29. 6612 = 3(q + 1) with q = 2203, a prime ≡ 3 (mod 8).

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I. This one is the case q = 2203 of the construction in
[PALEY_SKEW.md](../PALEY_SKEW.md), which gives a skew-Hadamard matrix of order 3(q + 1) for every prime power q ≡ 3 (mod
4).

- **Paley's matrix in sixteen blocks.** Paley's skew-Hadamard matrix of order q + 1 = 2204 on PG(1, 2203), signed and
  listed along the four orbits, of length s = 551, of g = (x ↦ −1/(x + 3))² as in
  [PALEY_BLOCKS.md](../../hadamard/PALEY_BLOCKS.md) (sign c = 1), is a 4 × 4 array of circulant s × s blocks N_kj.
- **The pairing.** After the base point of orbit 2 is moved by 173 (`shifts`), the blocks of column 0 in rows 1 and 2
  satisfy N_20 = N_10ᵀ. The pairing comes from an anti-automorphism of Paley's matrix that fixes two of the orbits and
  swaps the other two.
- **T-matrices of order 3.** T_0 = I, T_1 = Z, T_2 = −Z² and T_3 = 0, with Z the cyclic shift (`T`, one row per block
  row). The matrices A_j = Σ_k T_k ⊗ N_kj are ±1 matrices on ℤ₁₆₅₃ = ℤ₃ × ℤ₅₅₁, Σ_j A_j A_jᵀ = 6612 I, and A_0 = I ⊗
  N_00 + Z ⊗ N_10 − Z² ⊗ N_10ᵀ is skew-type.
- **The matrix.** The Goethals–Seidel array with A_0 first and R the back-diagonal identity is a skew-Hadamard matrix of
  order 6612. `sequences.txt` lists the four A_j as sequences on ℤ₁₆₅₃, A_0 first.

## Search

None: the matrix is determined by q. `python3 search/paley_blocks.py skew 2203` writes the certificate, and `python3
search/skew_3q.py 2203` builds a matrix of the same family from the closed formula of [PALEY_SKEW.md](../PALEY_SKEW.md)
and checks it.

## Status

**New** as a skew-Hadamard matrix, but obtainable classically ([PALEY_SKEW.md](../PALEY_SKEW.md#status)): Yamada's
sixteen-block form of Paley's matrix with its transposed pair (1991) and Seberry and Yamada's multiplication by
T-matrices (1990) were published, and we found them combined nowhere.

- The Đoković–Seberry census lists 1653 among its 737 undecided v, in none of its Paley, Szekeres, Szekeres–Whiteman or
  Spence series and not among its 41 sporadic values.
- Cati–Pasechnik Table 5 and Seberry–Yamada's Table 7.1 stop below v = 1000.

Doubling (Seberry) gives skew-Hadamard matrices of all orders 2ᵏ · 6612. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power and c (the blocks), shifts and columns (the pairing), the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences, one per line |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py skew/6612
```

SHA-256 of `H6612.csv`: `e2c652794da55bd2c1d8d9a35e49b22aa9d609ae75a34adcf21622ac2652e10a`

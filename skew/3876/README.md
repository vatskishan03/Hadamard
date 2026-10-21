# Skew-Hadamard matrix of order 3876

3876 = 4 · 969, with 969 = 3 · 17 · 19. 3876 = 3(q + 1) with q = 1291, a prime ≡ 3 (mod 8).

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I. This one is the case q = 1291 of the construction in
[PALEY_SKEW.md](../PALEY_SKEW.md), which gives a skew-Hadamard matrix of order 3(q + 1) for every prime power q ≡ 3 (mod
4).

- **Paley's matrix in sixteen blocks.** Paley's skew-Hadamard matrix of order q + 1 = 1292 on PG(1, 1291), signed and
  listed along the four orbits, of length s = 323, of g = (x ↦ −1/(x + 4)) as in
  [PALEY_BLOCKS.md](../../hadamard/PALEY_BLOCKS.md) (sign c = 1), is a 4 × 4 array of circulant s × s blocks N_kj.
- **The pairing.** After the base point of orbit 3 is moved by 171 (`shifts`), the blocks of column 0 in rows 2 and 3
  satisfy N_30 = N_20ᵀ. The pairing comes from an anti-automorphism of Paley's matrix that fixes two of the orbits and
  swaps the other two.
- **T-matrices of order 3.** T_0 = I, T_2 = Z, T_3 = −Z² and T_1 = 0, with Z the cyclic shift (`T`, one row per block
  row). The matrices A_j = Σ_k T_k ⊗ N_kj are ±1 matrices on ℤ₉₆₉ = ℤ₃ × ℤ₃₂₃, Σ_j A_j A_jᵀ = 3876 I, and A_0 = I ⊗ N_00
  + Z ⊗ N_20 − Z² ⊗ N_20ᵀ is skew-type.
- **The matrix.** The Goethals–Seidel array with A_0 first and R the back-diagonal identity is a skew-Hadamard matrix of
  order 3876. `sequences.txt` lists the four A_j as sequences on ℤ₉₆₉, A_0 first.

## Search

None: the matrix is determined by q. `python3 search/paley_blocks.py skew 1291` writes the certificate, and `python3
search/skew_3q.py 1291` builds a matrix of the same family from the closed formula of [PALEY_SKEW.md](../PALEY_SKEW.md)
and checks it.

## Status

**New, but obtainable classically** as a skew-Hadamard matrix ([PALEY_SKEW.md](../PALEY_SKEW.md#status)): Yamada's
sixteen-block form of Paley's matrix with its transposed pair (1991) and Seberry and Yamada's multiplication by
T-matrices (1990) were published, and we found them combined nowhere.

- The Đoković–Seberry census lists 969 among its 737 undecided v, in none of its Paley, Szekeres, Szekeres–Whiteman or
  Spence series and not among its 41 sporadic values.
- Cati–Pasechnik Table 5 gives 969(4): the smallest known skew order with odd part 969 was 16 · 969 = 15504.
  Seberry–Yamada's Table 7.1 (1992) agrees.

Doubling (Seberry) gives skew-Hadamard matrices of all orders 2ᵏ · 3876. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power and c (the blocks), shifts and columns (the pairing), the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences, one per line |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py skew/3876
```

SHA-256 of `H3876.csv`: `e9a8646a6df10a69e4602d030ff9cf7b73a3e532a68fa35fb4af2777d03477cc`

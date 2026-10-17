# Skew-Hadamard matrix of order 2076

2076 = 4 · 519, with 519 = 3 · 173. 2076 = 3(q + 1) with q = 691, a prime ≡ 3 (mod 8).

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I. This one is the case q = 691 of the construction in
[PALEY_SKEW.md](../PALEY_SKEW.md), which gives a skew-Hadamard matrix of order 3(q + 1) for every prime power q ≡ 3 (mod
4).

- **Paley's matrix in sixteen blocks.** Paley's skew-Hadamard matrix of order q + 1 = 692 on PG(1, 691), signed and
  listed along the four orbits, of length s = 173, of g = (x ↦ −1/(x + 4)) as in
  [PALEY_BLOCKS.md](../../hadamard/PALEY_BLOCKS.md) (sign c = 1), is a 4 × 4 array of circulant s × s blocks N_kj.
- **The pairing.** After the base point of orbit 2 is moved by 88 (`shifts`), the blocks of column 0 in rows 1 and 2
  satisfy N_20 = N_10ᵀ. The pairing comes from an anti-automorphism of Paley's matrix that fixes two of the orbits and
  swaps the other two.
- **T-matrices of order 3.** T_0 = I, T_1 = Z, T_2 = −Z² and T_3 = 0, with Z the cyclic shift (`T`, one row per block
  row). The matrices A_j = Σ_k T_k ⊗ N_kj are ±1 matrices on ℤ₅₁₉ = ℤ₃ × ℤ₁₇₃, Σ_j A_j A_jᵀ = 2076 I, and A_0 = I ⊗ N_00
  + Z ⊗ N_10 − Z² ⊗ N_10ᵀ is skew-type.
- **The matrix.** The Goethals–Seidel array with A_0 first and R the back-diagonal identity is a skew-Hadamard matrix of
  order 2076. `sequences.txt` lists the four A_j as sequences on ℤ₅₁₉, A_0 first.

## Search

None: the matrix is determined by q. `python3 search/paley_blocks.py skew 691` writes the certificate, and `python3
search/skew_3q.py 691` builds a matrix of the same family from the closed formula of [PALEY_SKEW.md](../PALEY_SKEW.md)
and checks it.

## Status

**New** as a skew-Hadamard matrix, but obtainable classically ([PALEY_SKEW.md](../PALEY_SKEW.md#status)): Yamada's
sixteen-block form of Paley's matrix with its transposed pair (1991) and Seberry and Yamada's multiplication by
T-matrices (1990) were published, and we found them combined nowhere.

- The Đoković–Seberry census lists 519 among its 737 undecided v, in none of its Paley, Szekeres, Szekeres–Whiteman or
  Spence series and not among its 41 sporadic values.
- Cati–Pasechnik Table 5 gives 519(4): the smallest known skew order with odd part 519 was 16 · 519 = 8304.
  Seberry–Yamada's Table 7.1 (1992) agrees.
- An ordinary Hadamard matrix of order 2076 is in [hadamard/2076](../../hadamard/2076) (Paley's matrix for the same q in
  sixteen blocks times T-matrices of order 3, without the pairing; First explicit construction).

Doubling (Seberry) gives skew-Hadamard matrices of all orders 2ᵏ · 2076. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power and c (the blocks), shifts and columns (the pairing), the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences, one per line |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py skew/2076
```

SHA-256 of `H2076.csv`: `5ec8970e369564c070acd357ec9f0d2b862b6ee091c4cc0fd4eed7e0577ef2dd`

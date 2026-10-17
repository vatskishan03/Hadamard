# Skew-Hadamard matrix of order 1692

1692 = 4 · 423, with 423 = 3² · 47. 1692 = 3(q + 1) with q = 563, a prime ≡ 3 (mod 8).

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I. This one is the case q = 563 of the construction in
[PALEY_SKEW.md](../PALEY_SKEW.md), which gives a skew-Hadamard matrix of order 3(q + 1) for every prime power q ≡ 3 (mod
4).

- **Paley's matrix in sixteen blocks.** Paley's skew-Hadamard matrix of order q + 1 = 564 on PG(1, 563), signed and
  listed along the four orbits, of length s = 141, of g = (x ↦ −1/(x + 6))² as in
  [PALEY_BLOCKS.md](../../hadamard/PALEY_BLOCKS.md) (sign c = 1), is a 4 × 4 array of circulant s × s blocks N_kj.
- **The pairing.** After the base point of orbit 3 is moved by 37 (`shifts`), the blocks of column 0 in rows 1 and 3
  satisfy N_30 = N_10ᵀ. The pairing comes from an anti-automorphism of Paley's matrix that fixes two of the orbits and
  swaps the other two.
- **T-matrices of order 3.** T_0 = I, T_1 = Z, T_3 = −Z² and T_2 = 0, with Z the cyclic shift (`T`, one row per block
  row). The matrices A_j = Σ_k T_k ⊗ N_kj are ±1 matrices on ℤ₃ × ℤ₁₄₁, Σ_j A_j A_jᵀ = 1692 I, and A_0 = I ⊗ N_00 + Z ⊗
  N_10 − Z² ⊗ N_10ᵀ is skew-type.
- **The matrix.** The Goethals–Seidel array with A_0 first and R: x ↦ −x is a skew-Hadamard matrix of order 1692.
  `sequences.txt` lists the four A_j as sequences on ℤ₃ × ℤ_s, at index i s + a, A_0 first.

## Search

None: the matrix is determined by q. `python3 search/paley_blocks.py skew 563` writes the certificate, and `python3
search/skew_3q.py 563` builds a matrix of the same family from the closed formula of [PALEY_SKEW.md](../PALEY_SKEW.md)
and checks it.

## Status

**New** as a skew-Hadamard matrix ([PALEY_SKEW.md](../PALEY_SKEW.md#status)).

- The Đoković–Seberry census lists 423 among its 737 undecided v, in none of its Paley, Szekeres, Szekeres–Whiteman or
  Spence series and not among its 41 sporadic values.
- Cati–Pasechnik Table 5 gives 423(4): the smallest known skew order with odd part 423 was 16 · 423 = 6768.
  Seberry–Yamada's Table 7.1 (1992) agrees, and Table 1.51 of the Handbook (2007) lists it as open.

Doubling (Seberry) gives skew-Hadamard matrices of all orders 2ᵏ · 1692. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power and c (the blocks), shifts and columns (the pairing), the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences, one per line |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py skew/1692
```

SHA-256 of `H1692.csv`: `8479c1177d9444811db403265fbfa793c0389d77c5191ead7b4127d65461d5e5`

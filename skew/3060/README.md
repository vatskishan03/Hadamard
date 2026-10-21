# Skew-Hadamard matrix of order 3060

3060 = 4 · 765, with 765 = 3² · 5 · 17. 3060 = 3(q + 1) with q = 1019, a prime ≡ 3 (mod 8).

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I. This one is the case q = 1019 of the construction in
[PALEY_SKEW.md](../PALEY_SKEW.md), which gives a skew-Hadamard matrix of order 3(q + 1) for every prime power q ≡ 3 (mod
4).

- **Paley's matrix in sixteen blocks.** Paley's skew-Hadamard matrix of order q + 1 = 1020 on PG(1, 1019), signed and
  listed along the four orbits, of length s = 255, of g = (x ↦ −1/(x + 6))² as in
  [PALEY_BLOCKS.md](../../hadamard/PALEY_BLOCKS.md) (sign c = 1), is a 4 × 4 array of circulant s × s blocks N_kj.
- **The pairing.** After the base point of orbit 2 is moved by 69 (`shifts`), the blocks of column 0 in rows 1 and 2
  satisfy N_20 = −N_10ᵀ. The pairing comes from an anti-automorphism of Paley's matrix that fixes two of the orbits and
  swaps the other two.
- **T-matrices of order 3.** T_0 = I, T_1 = Z, T_2 = Z² and T_3 = 0, with Z the cyclic shift (`T`, one row per block
  row). The matrices A_j = Σ_k T_k ⊗ N_kj are ±1 matrices on ℤ₃ × ℤ₂₅₅, Σ_j A_j A_jᵀ = 3060 I, and A_0 = I ⊗ N_00 + Z ⊗
  N_10 − Z² ⊗ N_10ᵀ is skew-type.
- **The matrix.** The Goethals–Seidel array with A_0 first and R: x ↦ −x is a skew-Hadamard matrix of order 3060.
  `sequences.txt` lists the four A_j as sequences on ℤ₃ × ℤ_s, at index i s + a, A_0 first.

## Search

None: the matrix is determined by q. `python3 search/paley_blocks.py skew 1019` writes the certificate, and `python3
search/skew_3q.py 1019` builds a matrix of the same family from the closed formula of [PALEY_SKEW.md](../PALEY_SKEW.md)
and checks it.

## Status

**New, but obtainable classically** as a skew-Hadamard matrix ([PALEY_SKEW.md](../PALEY_SKEW.md#status)): Yamada's
sixteen-block form of Paley's matrix with its transposed pair (1991) and Seberry and Yamada's multiplication by
T-matrices (1990) were published, and we found them combined nowhere.

- The Đoković–Seberry census lists 765 among its 737 undecided v, in none of its Paley, Szekeres, Szekeres–Whiteman or
  Spence series and not among its 41 sporadic values.
- Cati–Pasechnik Table 5 gives 765(4), and Seberry–Yamada's Table 7.1 (1992) agrees, but the family (|t| + 1)(q + 1)
  (Seberry Wallis) already gives 8 · 765 = 6120: q = 509 = 5² + 4 · 11², |t| + 1 = 12. So the smallest known skew
  order with odd part 765 was 8 · 765, and 4 · 765 was open.

Doubling (Seberry) gives skew-Hadamard matrices of all orders 2ᵏ · 3060. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, τ, power and c (the blocks), shifts and columns (the pairing), the T-matrices, the four sequences and the SHA-256 of the matrix |
| `sequences.txt` | the four sequences, one per line |

The matrix itself is not stored; `verify/verify.py` rebuilds it.

## Verification

```
python3 verify/verify.py skew/3060
```

SHA-256 of `H3060.csv`: `8e35bddc0c83212a67df1f7d7f5ef036c8995952772c192fb51258b60f0a4cc0`

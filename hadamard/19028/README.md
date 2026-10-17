# Hadamard matrix of order 19028

19028 = 4 · 4757, with 4757 = 67 · 71.

## Construction

Cooper–Wallis with Williamson-type matrices that are not Williamson matrices, as in [27604](../27604). The T-matrices
of order 71, in [`t-matrices/T71.txt`](../../t-matrices/T71.txt), are T-sequences: four sequences with entries 0, ±1
whose supports partition the positions and whose aperiodic autocorrelations sum to 0 off the origin. Their periodic
autocorrelations then sum to 0 as well, so the circulant matrices T₁, …, T₄ satisfy Σ TᵢTᵢᵀ = 71 I and give an
orthogonal design OD(284; 71, 71, 71, 71) through the Goethals–Seidel array. They are built by the classical route:

- Turyn-type sequences (X; Y; Z; W) of lengths 24, 24, 24, 23 with N_X + N_Y + 2N_Z + 2N_W = 0. These are taken from
  the table of Best, Đoković, Kharaghani and Ramp (2013), as stored in SageMath's `t_sequences` module
  (`turyn_type_sequences_smallcases(24)`).
- From them, base sequences (Z, W; Z, −W; X; Y) of lengths 47, 47, 24, 24.
- From those, T-sequences ((A + B)/2, 0₂₄), ((A − B)/2, 0₂₄), (0₄₇, (C + D)/2) and (0₄₇, (C − D)/2) of length 71.

Each step was checked exactly.

The Williamson-type matrices of order 67 come from the near-Williamson matrices in
[`near-williamson/NW67.txt`](../../near-williamson/NW67.txt) (Table 7 of arXiv:2605.08661v2, Kharaghani, Mohammadian,
Tayfeh-Rezaie and Zaitsev): circulant ±1 matrices A, B, C, D with B, C, D symmetric and ΣXXᵀ = 268 I, and
(AR, B, C, D) are Williamson-type. Substituting them for the variables of the OD gives a Hadamard matrix of order
4 · 71 · 67 = 19028. `build_cooper_wallis_near_williamson` in `verify/verify.py` checks both ingredients and
H Hᵀ = 19028 I.

## Status

**New, beyond the tables**, but obtainable classically (found by the novelty re-audit of October 16, after we had
labelled the order). Seberry–Yamada (1990, Theorem 4) applied to Paley's matrix of order 284 (q = 283 ≡ 3 mod 8), which
is an array of sixteen circulant blocks of order 71 (Seberry–Yamada 1993, Theorem 1), and to T-matrices of order 67
(Sawade 1985) gives the order from ingredients of 1993 or earlier ([PALEY_BLOCKS.md](../PALEY_BLOCKS.md)). We found no
construction or statement of the order itself.

The order lies above the range of Cati–Pasechnik Table 4, which stops at n = 2999, and is listed as a candidate in
[NOTES.md](../../NOTES.md) §7 ("4·4757 = 4·67·71 … using T-sequences of lengths 47 and 71"). The construction here is
that one.

- Paley I and II need 19027 = 53 · 359 and 9513 = 3² · 7 · 151 to be prime powers.
- 4757 is not a prime power (Miyamoto).
- 4755 = 3 · 5 · 317 is not a prime power (Yamada, Kiyasu).
- 4757 is not of Spence's form q² + q + 1.
- Cooper–Wallis needs T-sequences or T-matrices of length t and Williamson-type matrices of order w with
  tw = 4757.
  - t = 67, w = 71. Williamson-type matrices of order 71 are blank in Table A.1 of Seberry and Yamada (1992). The 2026
    near-Williamson paper establishes existence for every odd order up to 69, plus a list of larger orders that does
    not include 71.
  - t = 71, w = 67. This is the route used here. Williamson-type matrices of order 67 are listed in Table A.1 only as
    an 8-Williamson entry; four of order 67 appeared in 2026.
- We found no construction or statement of order 19028.

The ingredients are not new, only their combination. The sources searched are listed in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | T-matrices of order 71, near-Williamson matrices of order 67, SHA-256 |

The matrix is not stored.

## Verification

```
python3 verify/verify.py 19028
```

This rebuilds the matrix and checks it. It takes 54 seconds and about 1 GB of memory.

SHA-256 of `H19028.csv`: `ad9aec33fe30cb5de7e2d1c8fc2efaae04ec57c86c4d310362a2b28d4d2dbaac`

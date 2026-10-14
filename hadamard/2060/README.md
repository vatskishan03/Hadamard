# Hadamard matrix of order 2060

2060 = 4 · 103 · 5.

## Construction

Cooper–Wallis. Four T-matrices of order 103 are circulant matrices T₁, …, T₄ with entries 0, ±1. Their supports
partition ℤ₁₀₃, and Σ TᵢTᵢᵀ = 103 I. Putting them into the Goethals–Seidel array in the usual way gives an
orthogonal design OD(412; 103, 103, 103, 103). Substituting the four Williamson matrices of order 5 for its variables
gives a Hadamard matrix of order 2060.

`certificate.json` contains the first rows of the T-matrices (`0`, `+`, `-`) and of the Williamson matrices.

## Search

The T-matrices were found by a direct search:

- multiplier subgroup {1, 46, 56} of ℤ₁₀₃*,
- parity conditions solved over GF(2),
- exact meet-in-the-middle over the signs.

## Status

**Known.** Schneider posted a Hadamard matrix of order 2060 on 2026-08-23, and JD-Jones-ASES/Hadamard-T decoded
T-matrices of order 103 from it on 2026-08-29. The construction here was obtained independently. A Miyamoto
construction for 4 · 515 was claimed in 1991; Cati and Pasechnik could not verify it. The sources searched are listed
in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | T-matrices of order 103, Williamson matrices of order 5, SHA-256 |
| `H2060.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py 2060
gzip -dc hadamard/2060/H2060.csv.gz | ./check_matrix
```

SHA-256 of `H2060.csv`: `a89be4dfe3dafa847ec857ede44107947dc3dc879f33db8816196f68187341c8`

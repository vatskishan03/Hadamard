# Hadamard matrix of order 20712

20712 = 8 · 2589, with 2589 = 3 · 863 and 863 a prime ≡ 3 (mod 4).

## Construction

Miyamoto's quaternary array in blocks, times T-matrices of order 3
([QUATERNARY_MIYAMOTO_BLOCKS.md](../QUATERNARY_MIYAMOTO_BLOCKS.md)):

1. **K, a quaternary Hadamard matrix of order 862 = q − 1 with circulant blocks of order 431.** It is Paley's
   quaternary matrix of order (1723 + 1)/2 ([complex-hadamard/PALEY.md](../../complex-hadamard/PALEY.md)), written in
   the basis of two orbits of an element of order 431 of the non-split torus of SL(2, 1723)
   (`search/cmblocks/paley_circ.py 1723` reproduces the four first rows of the certificate).
2. **Sixteen bordered quaternary blocks of order q = 863** (points of GF(863)* along the primitive root 5), with the
   border matrices A and B of the certificate. Together they form a quaternary Hadamard matrix of order 4 · 863 =
   3452, an order with no known real Hadamard matrix, and they commute and are transposed by one permutation.
3. **Turyn's doubles of the blocks** (order 1726), which commute and are transposed by one signed permutation.
4. **T-matrices of order 3** (I, S, S², 0, S the cyclic shift) in the Goethals–Seidel array: order 8 · 863 · 3 = 20712.

`verify/verify.py` rebuilds K from its four first rows and checks K K* = 862 I, checks the border identities and the
T-matrices, assembles the matrix and checks H Hᵀ = 20712 I.

## Search

None for the matrix itself. The theorem and the border data were found in round 52
([search/NOTES.md](../../search/NOTES.md) §105): a search over 4 × 4 unit matrices gave the borders, and the plug-in was
found after the linear (quaternary) Goethals–Seidel plug-in had failed.

| Code | Command | Randomness and caps | Result |
|---|---|---|---|
| `search/cmblocks/make_cert.py`, `verify/verify.py` construction `quaternary-miyamoto-blocks` | `cd search/cmblocks; python3 make_cert.py 1723 3 cert.json` | none, deterministic | `certificate.json` (without its SHA-256), matrix SHA-256 `b5450dd9…` |

## Status

Cati–Pasechnik's Table 4 records 2589(4): order 16 · 2589 known, 8 · 2589 not; 4 · 2589 = 10356 and 4 · 863 = 3452
are open ([FRONTIER.md](../../FRONTIER.md)), so neither doubling nor Agayan's multiplication gives this order. The
literature check and an independent review are in progress (27 October); the label will be recorded in
[NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the first rows of the four circulant blocks of K (digits e for iᵉ), the borders A and B, the T-matrices, the SHA-256 of the matrix |

The matrix (20712² entries, about 860 MB as CSV) is not stored; `verify/verify.py` rebuilds it in a few minutes with
about 1.5 GB of memory. `search/cmblocks/rebuild.py` rebuilds it independently and also checks the intermediate claims.

## Verification

```
python3 verify/verify.py 20712
```

SHA-256 of `H20712.csv`: `b5450dd959630c2a9b337198e38d5a1b75b0239d5d078b8f216195f73ee3c61a`

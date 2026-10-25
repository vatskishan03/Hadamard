# Hadamard matrix of order 12952

12952 = 8 · 1619, with 1619 a prime ≡ 3 (mod 4).

## Construction

Seberry and Yamada's Hadamard matrices of order 8q (1992 survey, Lemma 8.12 and Corollaries 8.14, 8.18), with
symmetric Williamson-type matrices of order n = (q − 1)/2 = 809 that come from the symmetric Hadamard matrix of order
404 of this repository ([symmetric/404](../../symmetric/404)).

**Williamson-type matrices of order 809.** 809 is a prime ≡ 1 (mod 8) and (809 − 1)/2 = 404. Seberry and Yamada's
Corollary 8.8(ii) (Corollary 29(ii) of JCMCC 7, 1990) turns a symmetric Hadamard matrix of order 404 into symmetric
Williamson-type matrices of order 809, exactly as for [21036](../21036) (primitive root 3 mod 809).

**The 8q construction.** As for [8696](../8696), with q = 1619 and primitive root 2: eight pairwise amicable ±1 matrices
X_j of order q with Σ X_j X_jᵀ = 8q I, plugged into an OD(8; 1⁸). `verify/verify.py` rebuilds the symmetric matrix of
order 404, the Williamson-type matrices of order 809 and the X_j, checks each step exactly, and then H Hᵀ = 12952 I.

## Search

None beyond the symmetric matrix of order 404 ([symmetric/404](../../symmetric/404)). The route was found by the closure
computation of round 44 ([search/NOTES.md](../../search/NOTES.md) §95).

| Code | Command | Randomness and caps | Result |
|---|---|---|---|
| `search/sy8q.py`, `verify/verify.py` construction `seberry-yamada-8q` | `cd search; python3 sy8q.py c29 1619 symmetric/404/certificate.json` | none, deterministic | `certificate.json`, matrix SHA-256 `47d5f9c7…` |

## Status

**First explicit construction.** Cati–Pasechnik's Table 4 records 1619(3), that is order 8 · 1619 known. That entry
rests on Seberry–Yamada's Williamson-type matrices of order 809 (Table A.1, key wl), which come from Corollary 8.8(i),
whose printed proof does not give Williamson-type matrices ([NOTES.md](../../NOTES.md) §7). Part (ii) of the corollary,
used here, needs a symmetric Hadamard matrix of order 404 = 4 · 101, and 101 was undecided in the Đoković–Seberry
census until [symmetric/404](../../symmetric/404). We found no other valid construction:

- Paley's constructions do not apply: 12951 = 3² · 1439 and 6475 = 5² · 7 · 37 are not prime powers.
- No doubling: a Hadamard matrix of order 6476 = 4 · 1619 is not known (FRONTIER.md).
- Miyamoto's theorem for order 8q would need 2q − 3 = 3235 = 5 · 647 to be a prime power, and the complex Golay
  decomposition 8(x + y) would need two complex Golay lengths with sum 1619.

The order lies above the range of the 2026 survey's list of open orders (up to 10 000). The sources searched are listed
in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | q, g, the input certificate of order 404 with its SHA-256 and primitive root, the SHA-256 of the matrix |

The matrix (12952² entries) is not stored; `verify/verify.py` rebuilds it from the symmetric matrix of order 404.

## Verification

```
python3 verify/verify.py hadamard/12952
```

SHA-256 of `H12952.csv`: `47d5f9c72b35a271a39e455e3bb971b32618541f20a15c7600a3957392067f05`

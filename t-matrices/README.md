# T-matrices

T-matrices of order t are four circulant matrices with entries 0, ±1 whose supports
partition ℤ_t and with Σ TᵢTᵢᵀ = t I. By the Cooper–Wallis construction they give
an orthogonal design OD(4t; t, t, t, t). Replacing its variables with Williamson
matrices of order w gives a Hadamard matrix of order 4tw.

Each file `T<t>.txt` lists the first rows of T₁, …, T₄, one per line, written over
`0`, `+`, `-`. All of them were found by the search described in
[order 2060](../2060): a small multiplier subgroup, parity conditions over GF(2),
and meet-in-the-middle over the signs.

| t | Earlier appearances known to us |
|----:|---|
| 67 | classical (T-sequences of length 67); this instance from `attack/tmat_mult.py 67 29` (multiplier group of order 3) |
| 97 | none found; T-sequences of length 97 are open (Đoković, arXiv:1008.2043) |
| 103 | decoded from the 2060 matrix in JD-Jones-ASES/Hadamard-T (2026-08-29) |
| 109 | none found |
| 113 | Zuo–Xia–Xia (2013); Best–Đoković–Kharaghani–Ramp (2013) |
| 151 | posted by Schumer (2026-09-17) |
| 163 | JD-Jones-ASES/Hadamard-T (2026-08-29) |
| 181 | JD-Jones-ASES/Hadamard-T (2026-08-29) |
| 197 | none found |

For 97 and 109, which have the form 6m + 1, we could not check Lou, Xia and Xia,
"On T-matrices of order 6m + 1". Treat those two as possibly known.

To check a file, the conditions are that every column of the 4 × t array has exactly
one nonzero entry and that the periodic autocorrelations of the four rows sum to t
at shift 0 and to 0 elsewhere.

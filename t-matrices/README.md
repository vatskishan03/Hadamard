# T-matrices

T-matrices of order t are four circulant matrices with entries 0, ±1 whose supports
partition ℤ_t and with Σ TᵢTᵢᵀ = t I. By the Cooper–Wallis construction they give
an orthogonal design OD(4t; t, t, t, t). Replacing its variables with Williamson
matrices of order w gives a Hadamard matrix of order 4tw.

Each file `T<t>.txt` lists the first rows of T₁, …, T₄, one per line, written over
`0`, `+`, `-`. Orders 67, 101 and 131 come from `attack/tmat_mult.py`: a multiplier subgroup, a random split
of the orbits into the supports of T₁ + T₂ and T₃ + T₄, and an exact join of the four complete sign lists (numpy
for small lists, `attack/join4x.c` for large ones). The others were found by the search described in
[order 2060](../2060): a small multiplier subgroup, parity conditions over GF(2), and meet-in-the-middle over the
signs.

| t | Earlier appearances known to us |
|----:|---|
| 67 | classical (T-sequences of length 67); this instance from `attack/tmat_mult.py 67 29` (multiplier group of order 3) |
| 97 | none found; T-sequences of length 97 are open (Đoković, arXiv:1008.2043) |
| 101 | classical (T-sequences of length 1 + 10², Turyn 1974); this instance from `attack/tmat_mult.py 101 95` (multiplier group of order 5) |
| 103 | decoded from the 2060 matrix in JD-Jones-ASES/Hadamard-T (2026-08-29) |
| 109 | none found |
| 113 | Zuo–Xia–Xia (2013); Best–Đoković–Kharaghani–Ramp (2013) |
| 131 | none found; T-sequences of length 131 would need base sequences BS(66, 65) or Turyn-type sequences TT(44), neither known to us; from `attack/tmat_mult.py 131 53` (multiplier group of order 5, row sums 5, 0, −9, 5) |
| 151 | posted by Schumer (2026-09-17) |
| 163 | JD-Jones-ASES/Hadamard-T (2026-08-29) |
| 181 | JD-Jones-ASES/Hadamard-T (2026-08-29) |
| 197 | none found |

For 97 and 109, which have the form 6m + 1, we could not check Lou, Xia and Xia,
"On T-matrices of order 6m + 1". Treat those two as possibly known.

To check a file, the conditions are that every column of the 4 × t array has exactly
one nonzero entry and that the periodic autocorrelations of the four rows sum to t
at shift 0 and to 0 elsewhere.

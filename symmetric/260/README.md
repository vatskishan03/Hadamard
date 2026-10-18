# Symmetric Hadamard matrix of order 260

260 = 4 · 65, with 65 = 5 · 13: the smallest v of the Đoković–Seberry census that was undecided for symmetric matrices.

## Construction

A symmetric Hadamard matrix is a Hadamard matrix H with H = Hᵀ.

**The input.** H. Verma's quaternary Legendre pair (A, B) of length 64 (arXiv:2609.25573, 22 September 2026,
Table 1): sequences of fourth roots of unity on ℤ₆₄ with PAF_A + PAF_B = −2 at every nonzero shift. Two properties
that the paper does not mention make it usable here: A is real (its digits are 0 and 2, i.e. ±1), and B is symmetric,
B(−g) = B(g).

**Sequences.** Write B(g) = ((1 + i) p(g) + (1 − i) q(g))/2 with p, q = ±1 (digit 0: p = q = 1; 2: p = q = −1; 1:
p = 1, q = −1; 3: p = −1, q = 1). Since A is real, PAF_A + PAF_B = −2 becomes

    PAF_p + 2 PAF_A + PAF_q = −4   for every nonzero shift,

and p, q are symmetric because B is. So (p, A, A, q) is a bordered propus quadruple over ℤ₆₄ ([search/NOTES.md](../../search/NOTES.md)
§27), with row sums 2, 0, 0 (2² + 2 · 0² + 0² = 4). `sequences.txt` lists p, A and q.

**The matrix.** With the circulant matrices of p, A, q, the back-diagonal identity R, the propus array
P = [[p, A, A, q], [Aᵀ, qR, −p, −AR], [Aᵀ, −p, −qR, AR], [qᵀ, −AR, AR, −p]], K = H₂ ⊗ H₂ (Sylvester's Hadamard
matrix of order 4), S the same array on the row sums and L = −K S Kᵀ / 4,

    H = [[L, K ⊗ eᵀ], [Kᵀ ⊗ e, P]]

is a symmetric Hadamard matrix of order 4 · 64 + 4 = 260 (`verify.py` construction `propus-bordered`).

The quadruple has no symmetry beyond the reversal of p and q: no multiplier of ℤ₆₄ with a shift fixes all three
sequences up to sign. That is why the multiplier-invariant bordered propus searches of §27 and §28, which covered
every structure with at most 25–30 orbits on ℤ₆₄ and on the non-cyclic groups of order 64, could not find it.

## Search

None here: the sequences are a decoding of Verma's published pair (`python3 search/sym_qlp.py certify 65`). The
observation that the pair is of propus type was made by the exploration workflow of 17 October 2026; it was checked
from Verma's digit strings with three independent codes.

## Status

**New** as a symmetric Hadamard matrix. Ordinary Hadamard matrices of order 260 are classical; the symmetry is what
is new.

- The Đoković–Seberry census (2025-05-28) lists 65 among its 595 undecided v for symmetric matrices, as its smallest
  value. Đoković (arXiv:2101.05429v2) lists 65 first among the open v < 200, and the Handbook (2007, Table V.1.52, as
  quoted there) lists it as open.
- The input is from 2026. Verma states only the pair; Kotsireas and Winterhof's link from quaternary Legendre pairs
  gives quaternary and ordinary Hadamard matrices of orders 130 and 260, and no source we found notes that this pair
  gives a symmetric one.
- Zaitsev et al.'s quaternary matrix QH130 (page of 14 October 2026) has no anti-automorphism (nauty), so it is not
  equivalent to a symmetric matrix and does not give this order.
- Lebedev (arXiv:2609.04589) reports unsuccessful searches for quaternary Legendre pairs of length 64; other
  published pairs (lengths up to 34, Jedwab–Pender's families, Lebedev's lengths 42–80) give no remaining census value
  this way.

The sources searched are listed in [NOVELTY.md](../../NOVELTY.md).

## Files

| File | Contents |
|---|---|
| `certificate.json` | the three sequences p, A, q over ℤ₆₄, their row sums, K, and the SHA-256 of the matrix |
| `sequences.txt` | p, A and q, one per line |
| `H260.csv.gz` | the matrix |

## Verification

```
python3 verify/verify.py symmetric/260
gzip -dc symmetric/260/H260.csv.gz | ./check_matrix
```

`verify/verify.py` rebuilds the matrix from `certificate.json` (construction `propus-bordered`) and checks the symmetry
of p, the autocorrelation identity, H Hᵀ = 260 I and H = Hᵀ; `search/sym_qlp.py check` decodes Verma's pair and checks
PAF_A + PAF_B = −2, A real, B symmetric and the propus identity.

SHA-256 of `H260.csv`: `948586d6b68bc9a3c52ce2259903ad13f46ec54a800ee5cb8511b5a69d83eabe`

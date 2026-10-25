# tiedcore: exact tied joins over abelian groups

Programs of [../NOTES.md](../NOTES.md) §100. They search four ±1 sequences on a finite abelian group G that are
constant on the orbits of a multiplier group H, joined exactly on their periodic autocorrelations (PAF). G is a
product of finite fields F_q (additive group) and rings ℤ_k, written `F25xF29`, `Z2xF7xF17`, …; H acts factorwise
(a field factor by multiplication with a power of a primitive element, written as the exponent; a ring factor by a
unit). Python 3 with NumPy and SymPy.

| File | Contents |
|---|---|
| `gengroup.py` | `Group` (fields, rings, characters, multiplier action, subgroups of the multiplier group), `Engine` (exact integer PAF on the orbits of ⟨H, −1⟩ from the power spectrum at dual orbit representatives), `tied_run` (the zero-diagonal Szöllősi array, plain or bordered, with D = C ∘ n⁻¹), `full_check` (Σ PAF over the whole group by FFT) |
| `census.py` | the tied census of the complex Yamada–Kiyasu targets: every subgroup H with at most `maxC` C-orbits, A and B general or constant on quadratic or quartic classes, every class n of the multiplier group modulo H; resumable log |
| `gs_tied.py` | Goethals–Seidel quadruples (A, B, C, C ∘ n⁻¹) with A, B invariant under ⟨H, n⟩ |
| `gs_dtied.py` | Goethals–Seidel quadruples (X, X ∘ n⁻¹, Y, Y ∘ n⁻¹) |
| `anat850.py` | the anatomy of the order-850 core of [hadamard/6796](../../hadamard/6796) by field indices |

**The join.** For each class n, the list of PAF_A + PAF_B (A, B over their orbit signs) is matched against the list of
PAF_C + PAF_C ∘ n⁻¹ through a random linear key (int64, no overflow for the sizes used); every key match is re-checked
on the full integer PAF vector, and the first solution of each case is rebuilt over the whole group and checked by FFT.
The PAF of a sign vector is computed as |f̂|² at one character per dual orbit, mapped to the orbits of ⟨H, −1⟩ by the
exact inverse transform and rounded (the rounding error is asserted below 10⁻⁴).

**Checks.** The PAF matrix agrees with a direct FFT and with cyclotomic structure constants on random vectors.
Positive controls: GF(25) × GF(17) with H = ⟨(2, 6)⟩ gives the tied solutions of the order-850 core (512 for one
coarsening of A and B, 1536 over all classes n, 12 s); (5, 9) and (9, 13) give cores of orders 90 and 234; `gs_tied`
and `gs_dtied` find Goethals–Seidel quadruples over ℤ₃ × ℤ₅, F₃ × F₁₃ and F₅ × F₁₃.

```
python3 -c "from gengroup import *; G = Group('F25xF17'); print(tied_run(G, [(2, 6)], coarseA=[(4, 0), (0, 4)], coarseB=[(2, 0), (0, 2), (1, 1)])[0][0])"   # 512
python3 anat850.py                              # A, B, C, D of the order-850 core by (ind x, ind z)
python3 census.py 22 24 census.log F25xF29:0     # 261 cases, no solution (NOTES §100)
python3 gs_tied.py 20 24 gstied.log F5xF191     # 11 cases, no solution
```

Rules for long jobs: `nice`, `OMP_NUM_THREADS=1`, `ulimit -v 3500000`; the logs are resumable (one line per case).

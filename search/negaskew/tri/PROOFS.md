# The triangular family v = 1 + k(k + 1)/2: proofs of the structural facts of search/NOTES.md §92

Seven statements P1 to P7 about the negacyclic family of [NOTES.md](../../NOTES.md) §88 (skew-Hadamard matrices of order 8v from a Goethals–Seidel array with negacyclic blocks, v prime of the form 1 + k(k + 1)/2, k ≡ 3 (mod 4), blocks constant on the cosets of the subgroup of order k):
the parity of the autocorrelations, the row sums, the shape of the skew block, the line on which its spectrum lies, the symmetry group and the clusters of solutions, the period and Gauss-sum forms of the equations, and the absence of a congruence obstruction.
Each is PROVED here (P7 by an exact finite computation for the members listed) and used by the programs of this directory (the last section lists where). [`proofs_check.py`](proofs_check.py) tests every statement on the exact data (`python3 proofs_check.py`, 15 seconds, output at the end of this file); it checks, it does not prove.
P4(c) was first stated with the constants 2v − 5 = k² + k − 3 and 13 in place of 2v − 7 = k² + k − 5 and 21; the corrected arithmetic is below and the conclusion is unchanged (and now covers the exception v = 7 explicitly).

## Setting

v = 1 + k(k + 1)/2 is a prime, k ≡ 3 (mod 4), so e := (v − 1)/k = (k + 1)/2 is even; put h := e/2 (h is odd iff k ≡ 3 (mod 8)). ζ = e^{2πi/v}; g is a primitive root mod v, H = ⟨g^e⟩ is the subgroup of order k, C_j = g^j H (j ∈ ℤ_e) its cosets.
Since k is odd, −1 ∉ H; −1 = g^{kh} lies in C_h, so **−C_j = C_{j+h}**.

A *block* is a map u: ℤ_v → {1, i, −1, −i}, *H-invariant* here: u(0) = u₀ and u = c_j on C_j (a block is (u₀; c₀, …, c_{e−1})). Its *coset sum* is S = Σ_j c_j ∈ ℤ[i], its *row sum* R = Σ_x u(x) = u₀ + kS.
P_u(s) = Σ_x u(x + s)·conj u(x) is its periodic autocorrelation (PAF) and Û(t) = Σ_x u(x)ζ^{tx} its spectrum; P_u(−s) = conj P_u(s), P_u(as) = P_u(s) and |Û(at)| = |Û(t)| for a ∈ H, P_{λu} = P_u for |λ| = 1.
A *solution* is a quadruple of blocks A, B, C, D with Σ_X P_X(s) = 4v·δ_{s,0}, equivalently Σ_X |Û_X(t)|² = 4v for every t (the array of §88 is then Hadamard of order 8v); A is *twisted-skew*, u(−x) = i·conj u(x) for x ≠ 0 with u₀ ∈ {1, −i} (the array is skew), and B, C, D (the *free* blocks) have u₀ = 1 (a unit factor changes no PAF). The *normalised count* is the number of ordered solutions of this form (the COUNT of `qskew`).
Because of P_u(−s) = conj P_u(s) and the invariance under H, the PAF of a block is the *PAF vector* (Re P_u(g^t), Im P_u(g^t)) for t = 0, …, h − 1: e integers. Gaussian periods: η_j = Σ_{x∈C_j} ζ^x, so Σ_j η_j = −1 and η_{j+h} = conj η_j; the spectrum on the cosets is F_u(m) := Û(g^m) = u₀ + Σ_j c_j η_{j+m} (m ∈ ℤ_e).
Every unit c is written c = ((1 + i)p + (1 − i)q)/2 with p, q ∈ {±1} (1 = (+,+), i = (+,−), −1 = (−,−), −i = (−,+)): a block is a pair of ±1 sequences, as the ±1 negacyclic sequences of §88 are.

## P1. The parity of the PAF (PROVED)

*Statement.* Let n be odd and u: ℤ_n → {1, i, −1, −i} any map. For s ≢ 0 (mod n), Re P_u(s) is odd and Im P_u(s) is even.

*Proof.* Call x real if u(x) = ±1 and imaginary if u(x) = ±i. The term u(x + s)·conj u(x) is a unit that is real (±1) when x and x + s have the same type and imaginary (±i) when the types differ. The map x ↦ x + s splits ℤ_n into gcd(s, n) cycles of length n/gcd(s, n), and along a cycle the number of type changes is even.
So the number d of imaginary terms is even and the number n − d of real terms is odd. Re P_u(s) is a sum of n − d signs, hence odd; Im P_u(s) a sum of d signs, hence even. ∎

*Use.* `tri2` stops with "parity assumption violated" if a generated block breaks it, and takes (Re P − 1)/2 and (Im P)/2 as the integer coordinates of its hash; P7 uses it.

## P2. The row sums (PROVED for k ≥ 7)

*Statement.* Let k ≥ 7. In every solution exactly one block has S ≠ 0, and then S = u₀(1 ± i), R = u₀((k + 1) ± ki), |R|² = 4v − 3 = k² + (k + 1)²; the other three blocks have S = 0 and |R|² = 1. For the skew block S_A = (1 + i)m with m ≡ h (mod 2): **A is the big block iff k ≡ 3 (mod 8)** (then S_A = 1 + i with u₀ = 1, or S_A = −(1 + i) with u₀ = −i), and for k ≡ 7 (mod 8) A has S_A = 0 and one free block is big, with u₀ = 1 and S = 1 + i or 1 − i.

*Proof.* Summing Σ_X P_X(s) = 4vδ_{s,0} over all s gives Σ_X |R_X|² = 4v. Each coset sum S is a sum of e units; with S = a + bi the sum a + b has the parity of e, which is even, so S ∈ (1 + i)ℤ[i]: S = 0 or |S|² ≥ 2. If S ≠ 0 then |R| ≥ k|S| − 1 ≥ √2·k − 1; if S = 0 then |R| = |u₀| = 1.
Two blocks with S ≠ 0 would give Σ|R|² ≥ 2(√2 k − 1)² = 4k² − 4√2 k + 2 > 2k² + 2k + 4 = 4v for k ≥ 7; no block with S ≠ 0 gives Σ|R|² = 4 ≠ 4v. So exactly one block has S ≠ 0 and the other three have |R|² = 1, whence |u₀ + kS|² = 4v − 3 = 2k² + 2k + 1, that is 2 Re(conj(u₀)S) + k|S|² = 2k + 2.
Put |S|² = 2n (n a positive integer, since S ∈ (1 + i)ℤ[i]). Then Re(conj(u₀)S) = k + 1 − kn, and |Re(conj(u₀)S)| ≤ |S| = √(2n) forces n = 1 (for n ≥ 2, |Re(conj(u₀)S)| = k(n − 1) − 1 ≥ 7(n − 1) − 1 > √(2n)). So |S|² = 2 and Re(conj(u₀)S) = 1, i.e. conj(u₀)S = 1 ± i, S = u₀(1 ± i).
For the skew block, u(−x) = i·conj u(x) gives c_{j+h} = i·conj(c_j) for j < h, and c + i·conj c = ±(1 + i) for every unit c (1 ↦ 1 + i, i ↦ 1 + i, −1 ↦ −(1 + i), −i ↦ −(1 + i)). So S_A = (1 + i)m, m = Σ_{j<h} σ_j with σ_j = ±1, and m ≡ h (mod 2). For h odd (k ≡ 3 (mod 8)) m ≠ 0, so A carries the nonzero S: S_A = u₀(1 ± i) = (1 + i)m forces m = 1 and u₀ = 1, or m = −1 and u₀ = −i. For h even (k ≡ 7 (mod 8)) m is even, but a nonzero S has |S|² = 2, so m = 0 and S_A = 0; the big block is then free, with u₀ = 1. ∎

*Remarks.* (1) For k = 3 (v = 7) the statement fails: there are also solutions with |R|² = (13, 13, 1, 1) (`proofs_check.py P2`: 168 of the 248 ordered solutions at v = 7; the pattern (25, 1, 1, 1) has 80). (2) At k ≡ 3 (mod 8) the three free blocks have S = 0 and the number of such blocks is C(e, e/2)² (P3); at k ≡ 7 (mod 8) the big free block is counted by C(e, e/2 + 1)·C(e, e/2) for S = 1 + i.

## P3. The skew block is a pair of sign vectors (PROVED)

*Statement.* For a twisted-skew block put, for j < h, σ_j = (c_j + i·conj c_j)/(1 + i) and σ′_j = (c_j − i·conj c_j)/(1 − i). Both are ±1, the map c_j ↦ (σ_j, σ′_j) is the bijection 1 ↦ (+,+), i ↦ (+,−), −1 ↦ (−,−), −i ↦ (−,+), and c_{j+h} = i·conj(c_j). So a skew block is (u₀, σ, σ′) with σ, σ′ ∈ {±1}^h, and S_A = (1 + i)Σσ_j. With P2:
- for k ≡ 3 (mod 8): Σσ = 1 (u₀ = 1) or −1 (u₀ = −i), |A| = 2·C(h, (h + 1)/2)·2^h (8960 at v = 379, 640 at 191, 48 at 67);
- for k ≡ 7 (mod 8): Σσ = 0, u₀ ∈ {1, −i} free, |A| = 2·C(h, h/2)·2^h (2560 at v = 277, 16 at 29).
The free blocks with u₀ = 1 and S = 0 number |L| = C(e, e/2)² (11,778,624 at 379; 853,776 at 277; 63,504 at 191; 400 at 67), those with S = 1 + i number C(e, e/2 + 1)·C(e, e/2) (731,808 at 277).

*Proof.* The table c ↦ (c + i·conj c, c − i·conj c) = ((1 + i)σ, (1 − i)σ′) is a direct check on the four units, and it is invertible: c = ((1 + i)σ + (1 − i)σ′)/2. For c_{j+h} = i·conj c_j the pair (σ_j, σ′_j) of c_j determines the pair of c_{j+h}: it is (σ_j, −σ′_j). The sizes are counts of sign vectors with a prescribed sum times 2^h choices of σ′, times the two values of u₀ when both are allowed.
For a free block write c_j = ((1 + i)p_j + (1 − i)q_j)/2: its coset sum is S = ((1 + i)Σp + (1 − i)Σq)/2, so S = 0 iff Σp = Σq = 0 (C(e, e/2)² choices) and S = 1 + i iff Σp = 2, Σq = 0 (C(e, e/2 + 1)·C(e, e/2) choices). ∎

*Remark.* In the pair form the skew condition says p(−x) = p(x) and q(−x) = −q(x) for x ≠ 0 (p even, q odd), with q₀ = 1 and p₀ = 1 (u₀ = 1) or −1 (u₀ = −i); P4 and P6 read off the real and imaginary parts of the spectrum from this.

## P4. The spectrum of a skew block lies on a line; no flat skew block (PROVED; (c) corrected)

*Statement.* Let u be twisted-skew and t ∈ ℤ_v, including t = 0.
(a) conj Û(t) = (1 + i) − i·Û(t). Hence Û(t) = 1 + (1 + i)r_t with r_t real, |Û(t)|² = 2r_t² + 2r_t + 1, and ρ_t := (1 − i)Û(t) + i = 2r_t + 1 is real with ρ_t² = 2|Û(t)|² − 1.
(b) For t = g^m and a block with sign vectors σ, σ′: F(m) = u₀ + (1 + i)(x_m + y_m) with x_m = Σ_{j<h} σ_j X_{j+m}, y_m = Σ_{j<h} σ′_j Y_{j+m}, η_j = X_j + iY_j; ρ_m = 2(x_m + y_m) + ε with ε = 1 for u₀ = 1 and ε = −1 for u₀ = −i. Moreover Σ_{m∈ℤ_e} ρ_m = ε(e − 2) if the block is big (|R|² = 4v − 3) and ε·e if it is small (|R|² = 1), and Σ_m ρ_m² = e(2v − 7) in the first case, e(2v + 1) in the second.
(c) Let k ≥ 7. No t ≠ 0 and no twisted-skew u has |Û(t)|² = v − 3 or |Û(t)|² = v + 1. Call a block *flat* if |Û(t)|² is the same for all t ≠ 0; by Parseval the value is (v² − |R|²)/(v − 1), that is v + 1 for |R|² = 1 (PAF −1 off the origin; for v ≡ 3 (mod 4) the Legendre sequence with u₀ = 1 is such a block, and it is H-invariant) and v − 3 for |R|² = 4v − 3, the only row sums of P2. In a solution therefore (i) no twisted-skew block is flat, and (ii) the three free blocks are never all flat: for k ≡ 3 (mod 8) flat B, C, D give |Û_A(t)|² = 4v − 3(v + 1) = v − 3, for k ≡ 7 (mod 8) they give |Û_A(t)|² = 4v − (v − 3) − 2(v + 1) = v + 1.

*Proof.* (a) For x ≠ 0, u(−x) = i·conj u(x) gives conj u(x) = −i·u(−x). Then conj Û(t) = conj u₀ + Σ_{x≠0} conj u(x)ζ^{−tx} = conj u₀ − i Σ_{x≠0} u(−x)ζ^{−tx} = conj u₀ − i(Û(t) − u₀) = (conj u₀ + i u₀) − i·Û(t), and conj u₀ + i u₀ = 1 + i for both u₀ = 1 and u₀ = −i. Writing Û = a + bi, the real and imaginary parts of the identity both say a = 1 + b, so Û = 1 + b(1 + i), r_t = b. Then |Û|² = (1 + r)² + r² and ρ = (1 − i)(1 + (1 + i)r) + i = 1 + 2r.
(b) Û(g^m) = u₀ + Σ_j c_j η_{j+m}; with c_{j+h} = i·conj c_j, η_{j+h+m} = conj η_{j+m} and c_j = ((1 + i)σ_j + (1 − i)σ′_j)/2, i·conj c_j = ((1 + i)σ_j − (1 − i)σ′_j)/2, the pair (j, j + h) contributes (1 + i)σ_j X_{j+m} + (1 − i)σ′_j·iY_{j+m} = (1 + i)(σ_jX_{j+m} + σ′_jY_{j+m}). Then ρ_m = (1 − i)u₀ + 2(x_m + y_m) + i, which is 2(x_m + y_m) + 1 for u₀ = 1 and 2(x_m + y_m) − 1 for u₀ = −i.
Summing over m, Σ_m η_{j+m} = −1 gives Σ_m F(m) = e·u₀ − S, with S = (1 + i)m_A (P2): Σρ_m = (1 − i)(e u₀ − (1 + i)m_A) + ie = e − 2m_A for u₀ = 1 and −(e + 2m_A) for u₀ = −i, that is ε(e − 2|m_A|), with |m_A| = 1 (big) or 0 (small).
Parseval gives Σ_{t∈ℤ_v}|Û(t)|² = v², the t = 0 term is |R|² and the other t form e classes of k, so Σ_m |F(m)|² = (v² − |R|²)/k, that is e(v − 3) for |R|² = 4v − 3 and e(v + 1) for |R|² = 1 (v − 1 = ke). Then Σρ_m² = 2Σ|F(m)|² − e.
(c) Suppose |Û(t)|² = c ∈ {v − 3, v + 1} at some t ≠ 0. By (a) (2r_t + 1)² = N := 2c − 1, so N = 2v − 7 = k² + k − 5 or N = 2v + 1 = k² + k + 3, and 2r_t + 1 = ±√N. Now Û(t) lies in ℤ[i, ζ] ⊂ K = ℚ(i, ζ) = ℚ(ζ_{4v}), so √N is a real element of K.
If N is not a square, ℚ(√N) is a real quadratic subfield of K. The Galois group of ℚ(ζ_{4v}) is ℤ/2 × ℤ/(v − 1), which has three subgroups of index 2, so the quadratic subfields are ℚ(i), ℚ(√v*) and ℚ(√−v*) (v* = ±v ≡ 1 (mod 4)); exactly one of them is real, ℚ(√v). So N = v·a² with an integer a ≥ 1.
Case N = k² + k − 5. It is a square n² only if (2k + 1)² − 4n² = 21, i.e. (2k + 1 − 2n)(2k + 1 + 2n) = 21 = 1·21 = 3·7, i.e. k = 5 or k = 2. It equals v·a² only if v(2 − a²) = 7, i.e. a = 1 and v = 7, k = 3. Case N = k² + k + 3: n² gives 4n² − (2k + 1)² = 11, i.e. n = 3, k = 2; v·a² gives v(a² − 2) = 1, impossible. Neither case occurs for k ≡ 3 (mod 4), k ≥ 7. ∎

*Remarks.* (1) The exception is only v = 7, where 2v − 7 = 7 = v. (2) This excludes the ansatz B = C = D = Legendre (v ≡ 3 (mod 4), k ≡ 3 (mod 8)) for every member. A single flat free block is not excluded: it occurs at v = 29 and 67 (960 of 5952 and 576 of 6192 solutions, `stats2.py`) and its expected number is 61, 0.76 and 0.002 at 191, 379 and 631 (`legest.py`, HEURISTIC). (3) The three other spectra sum to 4v − |F_A(m)|² ≥ 0, so every A with a spectral peak max_m |F_A(m)|² > 4v is dead: 25.0%, 18.8% and 51.6% of the lists of A at v = 67, 191, 379 (`hyb.py`), to which the Gaussian estimate gives 5.7%, 5.4% and 16.6% of its mass (HEURISTIC; this contributes to, but cannot alone account for, the ratios 0.97, 0.49 and 0.28 of the exact or sampled count to the estimate; README.md).

## P5. The symmetries and the clusters of solutions (PROVED)

*Statement.* The following maps send solutions to solutions (normalised as above):
1. the multipliers x ↦ g^a x, a ∈ ℤ_e (ℤ_v*/H ≅ ℤ_e): all four blocks are shifted, c_j ↦ c_{j+a}; the reflection x ↦ −x is the shift by h;
2. the *tilde* of a free block, ũ(x) = conj u(−x) (c̃_j = conj c_{j+h}), independently for B, C, D; it keeps PAF, u₀ = 1 and H-invariance and sends S to conj S (so S = 1 + i to 1 − i);
3. the permutations of B, C, D;
4. **negA**: A′(x) = −A(x) for x ≠ 0 and A′(0) = −i·conj A(0) (u₀ = 1 ↔ −i); it is twisted-skew with S_{A′} = −S_A, and P_{A′} = P_A.
They generate the group Z_e × (Z₂ ≀ S₃) × Z₂ of order 96e (384, 576, 960, 1152, 1344 at e = 4, 6, 10, 12, 14). In a join with left pairs (A, B) and right pairs (C, D) of any linear residue of the PAF vector, the elements tilde_B, tilde_C, tilde_D, swap(C, D) and negA keep both pair sums; they generate a group of order 32 (order 16 when B is the big block with S = 1 + i, which tilde_B would send to S = 1 − i). So the solutions of one partition are a union of orbits of this group: **an orbit has 32/|Stab| elements (16/|Stab| in the second case); a solution without a stabiliser comes with 31 (15) others of the same left and right residue**, and the independent events of a sampling argument are the clusters, not the solutions.

*Proof.* Shifts: P_{u∘μ_a}(s) = Σ_x u(ax + as) conj u(ax) = P_u(as), so Σ_X P_{X∘μ_a} = (Σ_X P_X)∘μ_a = 0 off the origin; μ_a keeps u₀, the skew condition (u(−ax) = i·conj u(ax)) and the coset sums. Tilde: P_ũ(s) = Σ_x conj u(−x − s)·u(−x) = Σ_z u(z + s)·conj u(z) = P_u(s); the rest is read off the definition. Permutations: the equation is symmetric in B, C, D.
negA: write A = u₀δ₀ + a with a supported off the origin, A′ = u₀′δ₀ − a, u₀′ = −i·conj u₀. For s ≠ 0, P_A(s) = a(s)·conj u₀ + u₀·conj a(−s) + Σ_{x≠0,−s} a(x + s)·conj a(x), and P_{A′}(s) = −a(s)·conj u₀′ − u₀′·conj a(−s) + (the same sum). By a(−s) = i·conj a(s), conj a(−s) = −i·a(s), so P_A(s) − P_{A′}(s) = a(s)·[conj u₀ + conj u₀′ − i(u₀ + u₀′)] = 0 for (u₀, u₀′) = (1, −i) and for (−i, 1) (the bracket is (1 + i) − i(1 − i) = 0 in both cases). A′ is twisted-skew: A′(−x) = −a(−x) = −i·conj a(x) = i·conj A′(x).
The group generated has the stated order as an abstract group (tilde and the permutations generate Z₂ ≀ S₃ of order 48; shifts and negA commute with everything). The pair sums: PAF(tilde X) = PAF(X), PAF(negA) = PAF(A), and the swap exchanges the two terms of the right sum. ∎

*Consequences.* (1) At v = 379 (k ≡ 3 (mod 8), group of order 32) every hit partition of the sample holds exactly 32 solutions (14 of 14: one orbit each); at v = 277 (k ≡ 7 (mod 8), group of order 16) the hits per partition of the unreduced run are multiples of 8, 16 for 2322 of the 2996 hit partitions. (2) The classes (orbits of the full group of order 96e) at v = 29, 67, 191 number 24, 13, 105, at 277 313 (`classes_red.py`: 311 of size 1152 and 2 of size 576). At 191 the orbit sizes are 960 = 96e except for 3 orbits of size 480, and the only stabilisers found (`stab.py`) are single tildes, that is solutions with a self-dual free block (c_{j+h} = conj c_j); at 29 and 67 (e ≤ 6) there are also quarter orbits, from equal free blocks. (3) A run of 14 even residues of a join modulo 2^18 at v = 379 meets E/(32·2^17)·14 clusters, not E/2^18·14 hits (README.md).

## P6. The period form and the Gauss-sum form of the equations (PROVED)

*Statement.* The PAF equations of a solution are equivalent to Σ_X |R_X|² = 4v together with Σ_X |F_X(m)|² = 4v for m ∈ ℤ_e (one linear relation links these e + 1 equations, so e of them are independent). Put d_j = c_j + c_{j+h}, d′_j = c_j − c_{j+h} (j < h) and X_j = Re η_j, Y_j = Im η_j (X_{j+h} = X_j, Y_{j+h} = −Y_j).
*Period form.* For m ∈ ℤ_h, F(m) = U_m + iV_m and F(m + h) = U_m − iV_m with U_m = u₀ + Σ_{j<h} d_jX_{j+m} and V_m = Σ_{j<h} d′_jY_{j+m} (complex numbers); the equations are Σ_X (|U_X|² + |V_X|²) = 4v and Σ_X Im(U_X·conj V_X) = 0 for all m ∈ ℤ_h. For the skew block U = u₀ + (1 + i)x, V = (1 − i)y, with x, y of P4(b), and Im(U·conj V) = y(2x + 1) when u₀ = 1.
*Gauss-sum form.* Let ω = e^{2πi/e}, χ the character of ℤ_v* with χ(g) = ω, G_l = Σ_{x≠0} χ^l(x)ζ^x (l ∈ ℤ_e; G₀ = −1, |G_l|² = v for l ≠ 0), ĉ(l) = (1/e)Σ_j c_jω^{−lj}, B(l) = ĉ(l)G_l for l ≠ 0 and B(0) = u₀ − S/e. Then F(m) = Σ_l B(l)ω^{−lm} and the equations read Σ_X PAF^{ℤ_e}_{B_X}(s) = 4v·δ_{s,0} with PAF^{ℤ_e}_B(s) = Σ_l B(l + s)·conj B(l): **four sequences of length e whose periodic autocorrelations sum to 4v δ.**

*Proof.* The sum Σ_X P_X is H-invariant and its Fourier transform at t is Σ_X |Û_X(t)|²; it equals 4vδ₀ iff the transform is 4v for every t, and by H-invariance it suffices to take t = 0 and t = g^m. For any four blocks Parseval gives Σ_X |R_X|² + k Σ_m Σ_X |F_X(m)|² = 4v², which is the relation: once Σ_X |R_X|² = 4v holds, the e left sides at t = g^m add up to 4v·e, the sum of the right sides.
Period form: F(m) = u₀ + Σ_{j<h}[c_jη_{j+m} + c_{j+h}η_{j+h+m}] = u₀ + Σ_{j<h}[(c_j + c_{j+h})X_{j+m} + i(c_j − c_{j+h})Y_{j+m}] = U_m + iV_m, and with Y_{j+h+m} = −Y_{j+m}, F(m + h) = U_m − iV_m. Then |F(m)|² + |F(m + h)|² = 2(|U|² + |V|²) and |F(m)|² − |F(m + h)|² = 4 Im(U·conj V). For the skew block d_j = (1 + i)σ_j, d′_j = (1 − i)σ′_j (P3), so U = u₀ + (1 + i)x, V = (1 − i)y, and for u₀ = 1, Im(U·conj V) = Im((1 + (1 + i)x)(1 + i)y) = y + 2xy.
Gauss-sum form: the indicator of C_j on x ≠ 0 is (1/e)Σ_l χ^l(x)ω^{−lj}, so η_j = (1/e)Σ_l ω^{−lj}G_l and F(m) = u₀ + Σ_j c_jη_{j+m} = u₀ + Σ_l ω^{−lm}G_lĉ(l); the l = 0 term is −S/e. The ℤ_e-Fourier transform of |F(m)|² = Σ_{l,l′}B(l)conj B(l′)ω^{−(l−l′)m} is e·PAF^{ℤ_e}_B(s), and 4v is the transform of 4v·e δ. |G_l|² = v for a nontrivial character is the classical Gauss-sum identity. ∎

## P7. No congruence obstruction (EXACT for the members listed)

*Statement.* Let w(u) ∈ ℤ^{e−1} be the PAF vector of a block without the coordinate Re P(g^{h−1}) (given the row sums, that coordinate of Σ_X P_X is determined: Σ_{s≠0} P(s) = |R|² − v). For v = 29, 67, 191, 277, 379, 631, 947, 1129 the lattice Λ generated by the differences w(u) − w(u′) of blocks of the same type (skew A, small free, and for k ≡ 7 (mod 8) big free) has rank e − 1 and index 2^{e−1} in ℤ^{e−1}, so Λ = (2ℤ)^{e−1}. Consequently the coset of Λ that contains the sum Σ_X w(u_X) of four PAF vectors is Λ itself, which contains 0: no congruence on the PAF vector (modulo (1 + i)^r, 2^r or any other modulus) excludes a solution. (At v = 379: the lattice has Euclidean covolume 2^{13}√7 = 21,674 in the e coordinates, the number `tri_est.py` prints; for 631 it is 2^{17}·3 = 393,216.)

*Proof.* By P1 every coordinate of w(u) − w(u′) is even (real parts odd − odd, imaginary parts even − even), so Λ ⊆ (2ℤ)^{e−1}, which has index 2^{e−1}. The lattice generated by sampled differences of actual blocks (`proofs_check.py P7`: 800 samples per type, exact integer elimination) is contained in Λ and has index 2^{e−1}; hence it equals (2ℤ)^{e−1} and so does Λ. Every Σ_X w(u_X) lies in the coset s₀ + Λ of the sum s₀ = Σ_X w(u_X⁰) of reference blocks of the four types, and s₀ ∈ (2ℤ)^{e−1} = Λ (four vectors with odd real and even imaginary entries add up to even entries), so the coset is Λ itself. ∎

*Remark.* This is the reason the Gaussian estimate of `est_fast.py` uses the covolume 2^{e−1} without the lattice step; it excludes a parity or modular obstruction only, not a non-linear one.

## What is not proved

The existence of solutions for a given v (EXACT for 7, 29, 67, 191, 277, 379 by the joins, UNKNOWN for 631, 947, 1129, 1327, …), the exact count at 379, the Gaussian-estimate ratios and the finiteness of the family (HEURISTIC): README.md.

## Where the programs use the statements

| Statement | Program |
|---|---|
| P1 | `tri2.c` (`fill`: aborts if a block breaks the parity, takes (Re P − 1)/2 and Im P/2), `proofs_check.py` |
| P2 | `tri2.c` modes 3 and 7 (which block is big), `tri_join.c`, `ana.py` (the row-sum pattern of a solution file) |
| P3 | `tri2.c` (`genA`, `genSmall`, `genBig`), `astruct.py`, the list sizes in README.md |
| P4 | `leg.py`, `legest.py`, `stats2.py`, `flat.py`, `astat.py`, `hyb.py`, `cop.py` (the spectrum of A, the positivity cut-off) |
| P5 | `tri2.c` (`red=1`: the weights of the reduced hits), `symm.py`, `symm2.py`, `stab.py`, `classes_red.py`, `resid379.py` |
| P6 | `tricore.py` (spectra from the periods), `gauss_dom.py`, `est_fast.py` (the PAF coordinates) |
| P7 | `coset.py`, `est_all.py`, `est_fast.py`, `tri_est.py` |

## Output of `python3 proofs_check.py` (15 s wall; run in the staging tree after the builds of README.md)

```
P1  parity of the PAF
  ok   P1: Re P(s) odd and Im P(s) even for every s != 0   [70200 (sequence, shift) pairs over v = 7, 15, 21, 29, 67, 101]
P2  row sums
  ok   P2 at v = 29 (k = 7, k mod 8 = 7): every solution has |R|^2 = (4v-3, 1, 1, 1) in some order   [992 solutions, patterns {(1, 113, 1, 1): 992}]
  ok      A is small and one free block is big (k = 7 mod 8)
  ok   P2 at v = 67 (k = 11, k mod 8 = 3): every solution has |R|^2 = (4v-3, 1, 1, 1) in some order   [6192 solutions, patterns {(265, 1, 1, 1): 6192}]
  ok      A is the big block (k = 3 mod 8)
  ok   P2 on the recorded data/s379.out   [448 solutions, patterns {(1513, 1, 1, 1): 448}]
  ok   P2 on the recorded data/c277_red.out   [646 solutions, patterns {(1, 1105, 1, 1): 646}]
  ok   P2 needs k >= 7: at v = 7 (k = 3) the patterns are (25,1,1,1) and (13,13,1,1)   [{(25, 1, 1, 1): 80, (13, 13, 1, 1): 168}]
P3  the skew blocks
  ok   P3: c -> (sigma, sigma') is a bijection: 1 -> (+,+), i -> (+,-), -1 -> (-,-), -i -> (-,+)   [{0: (1, 1), 1: (1, -1), 2: (-1, -1), 3: (-1, 1)}]
  ok   P3: c = ((1+i) sigma + (1-i) sigma') / 2
  ok   P3 at v = 29 (h = 2): |A| = 16 by enumeration; |L| = C(e,e/2)^2 = 36, |Big| = C(e,e/2+1) C(e,e/2) = 24   [formula 16]
  ok   P3 at v = 67 (h = 3): |A| = 48 by enumeration; |L| = C(e,e/2)^2 = 400, |Big| = C(e,e/2+1) C(e,e/2) = 300   [formula 48]
  ok   P3 at v = 191 (h = 5): |A| = 640 by enumeration; |L| = C(e,e/2)^2 = 63504, |Big| = C(e,e/2+1) C(e,e/2) = 52920   [formula 640]
  ok   P3 at v = 277 (h = 6): |A| = 2560 by enumeration; |L| = C(e,e/2)^2 = 853776, |Big| = C(e,e/2+1) C(e,e/2) = 731808   [formula 2560]
  ok   P3 at v = 379 (h = 7): |A| = 8960 by enumeration; |L| = C(e,e/2)^2 = 11778624, |Big| = C(e,e/2+1) C(e,e/2) = 10306296   [formula 8960]
  ok   P3 at v = 631 (h = 9): |A| = 129024 by enumeration; |L| = C(e,e/2)^2 = 2363904400, |Big| = C(e,e/2+1) C(e,e/2) = 2127513960   [formula 129024]
  ok   P3: every A of the exact set of v = 29 satisfies c_{j+h} = i conj(c_j), u0 in {1, -i}
  ok   P3: every A of the exact set of v = 67 satisfies c_{j+h} = i conj(c_j), u0 in {1, -i}
P4  the spectrum of a skew block
  ok   P4(a) at v = 29: conj U(t) = (1+i) - iU(t) at every t (max error 2.4e-13); U = 1 + (1+i) r with r real; |U|^2 = 2r^2 + 2r + 1   [16 distinct A]
  ok   P4(b) at v = 29: rho_m = (1-i) F(m) + i real, sum rho_m = +- e, sum rho_m^2 = e(2v+1)   [e = 4]
  ok   P4(a) at v = 67: conj U(t) = (1+i) - iU(t) at every t (max error 8.7e-13); U = 1 + (1+i) r with r real; |U|^2 = 2r^2 + 2r + 1   [24 distinct A]
  ok   P4(b) at v = 67: rho_m = (1-i) F(m) + i real, sum rho_m = +-( e - 2), sum rho_m^2 = e(2v-7)   [e = 6]
  ok   P4(c): for k = 3 mod 4, 7 <= k <= 10^6: 2v-7 = k^2+k-5 and 2v+1 are neither a square nor v times a square   [[]]
  ok   P4(c): k^2 + k - 5 is a square only for k = 2, 5 (identity (2k+1)^2 - 4n^2 = 21); at v = 7 (k = 3) 2v - 7 = 7 = v * 1^2 is the exception   [squares at k = [2, 5]]
P5  symmetries and clusters
  ok   P5 at v = 29: the 15360 images of 40 solutions under the 384 elements of Z_e x tilde^3 x S_3 x negA are solutions
  ok      (mode 7: the images with the big block at B and S = 1+i, 2560 of them, lie in the list)
  ok      the left sum PAF(A)+PAF(B) and the right sum PAF(C)+PAF(D) are invariant under tilde_B, tilde_C, tilde_D, swap(C,D), negA (exact integer vectors)
  ok   P5 at v = 67: the 23040 images of 40 solutions under the 576 elements of Z_e x tilde^3 x S_3 x negA are solutions
  ok      (mode 3: the exact list is the whole normalised set, so every image lies in it)   [23040 of 23040]
  ok      the left sum PAF(A)+PAF(B) and the right sum PAF(C)+PAF(D) are invariant under tilde_B, tilde_C, tilde_D, swap(C,D), negA (exact integer vectors)
  ok   P5: each of the 14 hit partitions of the 379 sample is exactly one orbit of the group of order 32 (tilde B, C, D; swap C, D; A -> -A)   [448 solutions, orbit sizes [32]]
P6  period form and Gauss-sum form
  ok   P6 at v = 29: eta_{j+h} = conj eta_j, sum_j eta_j = -1, and F(m) = u0 + sum_j c_j eta_{j+m} equals the direct spectrum
  ok      sum_X |F_X(m)|^2 = 4v for every m   [max error 7.1e-14]
  ok      F(m) = U_m + iV_m, F(m+h) = U_m - iV_m; sum_X (|U|^2+|V|^2) = 4v and sum_X Im(U conj V) = 0 for m in Z_h   [max errors 2.8e-14, 8.0e-15]
  ok      skew block A: U_m = u0 + (1+i) x_m, V_m = (1-i) y_m, Im(U conj V) = y(2x+1) for u0 = 1
  ok      Gauss-sum form: G_0 = -1, |G_l|^2 = v, F(m) = sum_l B(l) w^(-lm), sum_X PAF^{Z_e}_{B_X}(s) = 4v delta_s   [max error 4.3e-14]
  ok   P6 at v = 67: eta_{j+h} = conj eta_j, sum_j eta_j = -1, and F(m) = u0 + sum_j c_j eta_{j+m} equals the direct spectrum
  ok      sum_X |F_X(m)|^2 = 4v for every m   [max error 4.5e-13]
  ok      F(m) = U_m + iV_m, F(m+h) = U_m - iV_m; sum_X (|U|^2+|V|^2) = 4v and sum_X Im(U conj V) = 0 for m in Z_h   [max errors 5.7e-14, 4.3e-14]
  ok      skew block A: U_m = u0 + (1+i) x_m, V_m = (1-i) y_m, Im(U conj V) = y(2x+1) for u0 = 1
  ok      Gauss-sum form: G_0 = -1, |G_l|^2 = v, F(m) = sum_l B(l) w^(-lm), sum_X PAF^{Z_e}_{B_X}(s) = 4v delta_s   [max error 2.8e-13]
  ok   P6 at v = 379: eta_{j+h} = conj eta_j, sum_j eta_j = -1, and F(m) = u0 + sum_j c_j eta_{j+m} equals the direct spectrum
  ok      sum_X |F_X(m)|^2 = 4v for every m   [max error 3.4e-11]
  ok      F(m) = U_m + iV_m, F(m+h) = U_m - iV_m; sum_X (|U|^2+|V|^2) = 4v and sum_X Im(U conj V) = 0 for m in Z_h   [max errors 9.1e-13, 3.8e-13]
  ok      skew block A: U_m = u0 + (1+i) x_m, V_m = (1-i) y_m, Im(U conj V) = y(2x+1) for u0 = 1
  ok      Gauss-sum form: G_0 = -1, |G_l|^2 = v, F(m) = sum_l B(l) w^(-lm), sum_X PAF^{Z_e}_{B_X}(s) = 4v delta_s   [max error 3.3e-11]
P7  the lattice of the PAF vectors
  ok   P7 at v = 29 (e = 4): differences even; rank 3, index 8 = 2^3 in Z^3
  ok   P7 at v = 67 (e = 6): differences even; rank 5, index 32 = 2^5 in Z^5
  ok   P7 at v = 191 (e = 10): differences even; rank 9, index 512 = 2^9 in Z^9
  ok   P7 at v = 277 (e = 12): differences even; rank 11, index 2048 = 2^11 in Z^11
  ok   P7 at v = 379 (e = 14): differences even; rank 13, index 8192 = 2^13 in Z^13
  ok   P7 at v = 631 (e = 18): differences even; rank 17, index 131072 = 2^17 in Z^17
  ok   P7 at v = 947 (e = 22): differences even; rank 21, index 2097152 = 2^21 in Z^21
  ok   P7 at v = 1129 (e = 24): differences even; rank 23, index 8388608 = 2^23 in Z^23
ALL CHECKS PASSED
```
